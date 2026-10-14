// Interactive ice-sheet map (OpenLayers): raster tiles in the polar
// stereographic projection of the posters, with vector overlays, name search
// and a legend. All data come from assets/<region>/ (scripts/web.jl).
(function () {
  const el = document.getElementById('map');
  const region = el.dataset.region;
  const base = 'assets/' + region + '/';
  fetch(base + 'config.json').then(r => r.json()).then(init);

  // label styles per type (after LSTYLE in scripts/labels.jl)
  const LSTYLE = {
    region:     {font: 'bold 14px sans-serif',   color: '#4d4d4d', upper: true},
    mountains:  {font: 'italic 13px sans-serif', color: '#8b4726'},
    basin:      {font: 'italic 12px sans-serif', color: '#595959'},
    dome:       {font: '12px sans-serif', color: '#000', marker: 'triangle'},
    station:    {font: '12px sans-serif', color: '#000', marker: 'circle'},
    icecore:    {font: '12px sans-serif', color: '#8b0000', marker: 'diamond'},
    lake:       {font: 'italic 12px sans-serif', color: '#000080', marker: 'circle'},
    island:     {font: '12px sans-serif', color: '#333'},
    sea:        {font: 'italic 16px sans-serif', color: '#1d4e89', halo: false},
    settlement: {font: '12px sans-serif', color: '#000', marker: 'circle'},
    glacier:    {font: '12px sans-serif', color: '#000', marker: 'dot'},
    iceshelf:   {font: 'bold 12px sans-serif', color: '#000'},
    fjord:      {font: 'italic 12px sans-serif', color: '#000080'},
  };
  // finest resolution (m/px) at which tier 2 names and gazetteer names appear
  const TIER2_RES = 1500, GAZ_RES = 700;

  function marker(kind, color) {
    const fill = new ol.style.Fill({color}), stroke = new ol.style.Stroke({color: '#fff', width: 1.2});
    if (kind === 'triangle') return new ol.style.RegularShape({points: 3, radius: 6, fill, stroke});
    if (kind === 'diamond') return new ol.style.RegularShape({points: 4, radius: 6, fill, stroke});
    if (kind === 'dot') return new ol.style.Circle({radius: 2.5, fill: new ol.style.Fill({color: '#000'})});
    return new ol.style.Circle({radius: 4.5, fill, stroke});
  }

  const styleCache = {};
  function labelStyle(f) {
    const type = f.get('type'), name = f.get('name');
    const key = type + '|' + name;
    if (styleCache[key]) return styleCache[key];
    const s = LSTYLE[type] || LSTYLE.glacier;
    const text = new ol.style.Text({
      text: s.upper ? name.toUpperCase() : name, font: s.font,
      fill: new ol.style.Fill({color: s.color}),
      stroke: s.halo === false ? undefined : new ol.style.Stroke({color: 'rgba(255,255,255,0.85)', width: 3}),
      textAlign: s.marker ? 'left' : 'center', offsetX: s.marker ? 8 : 0,
    });
    return (styleCache[key] = new ol.style.Style({image: s.marker ? marker(s.marker, s.color) : undefined, text}));
  }

  function lineStyle(color, width, lineDash) {
    return new ol.style.Style({stroke: new ol.style.Stroke({color, width, lineDash})});
  }

  function fmt(v, d) { return v.toLocaleString('en', {maximumFractionDigits: d, minimumFractionDigits: d}); }

  // key numbers of the whole ice sheet and, where given (Antarctica), of its regions
  function numbersTable(n) {
    const cols = [['All', '']].concat(n.east_ice_area_km2 === undefined ? [] :
      [['East', 'east_'], ['West', 'west_'], ['Pen.', 'peninsula_']]);
    const rows = [['Ice area [10⁶ km²]', 'ice_area_km2', 1e6, 2], ['&nbsp;floating [10³ km²]', 'floating_area_km2', 1e3, 0],
                  ['Ice volume [10⁶ km³]', 'ice_volume_km3', 1e6, 2], ['Max. thickness [m]', 'max_thickness_m', 1, 0],
                  ['Sea-level equiv. [m]', 'sea_level_equivalent_m', 1, 1]];
    const head = cols.length > 1 ? `<tr><th></th>${cols.map(c => `<th>${c[0]}</th>`).join('')}</tr>` : '';
    return `<table class="cryo-numbers">${head}${rows.map(([name, key, f, d]) =>
      `<tr><td>${name}</td>${cols.map(c => `<td>${fmt(n[c[1] + key] / f, d)}</td>`).join('')}</tr>`).join('')}</table>` +
      (cols.length > 1 ? '<small>Regions of IMBIE 2; ice shelves and islands go to the nearest region.</small>' : '');
  }


  // --- adaptive colour scales ------------------------------------------------
  // The surface and bed styles are tiles of packed data (R: hillshade, G/B: class and elevation,
  // scripts/web.jl `encode_raster`), coloured here with limits taken from the part of the map in
  // view. The colour maps are those of the posters (`compose` in scripts/common.jl).
  const TILE_PX = 256, HIST_N = 16384, ADAPT_MS = 120, RAW_CACHE = 256;
  const clamp01 = t => t < 0 ? 0 : t > 1 ? 1 : t;
  const hexRgb = h => [1, 3, 5].map(k => parseInt(h.slice(k, k + 2), 16));

  // lookup tables [r, g, b, ...] of the 256-step colour maps
  function rampTables(ramps) {
    const out = {};
    for (const [k, stops] of Object.entries(ramps)) {
      out[k] = new Uint8Array(3 * stops.length);
      stops.forEach((h, i) => out[k].set(hexRgb(h), 3 * i));
    }
    return out;
  }

  // colour of one packed pixel into out[o..o+3] (after `compose`: colour map, hillshade, then mixing)
  function paintTile(raw, mode, lim, R, T, out) {
    const shelf = R.shelfRgb, n = T.ice.length / 3 - 1;
    const srf = lim.srf, ocean = lim.ocean, bed = lim.bed;
    for (let o = 0; o < raw.length; o += 4) {
      const r = raw[o], g = raw[o + 1], b = raw[o + 2];
      let c0, c1, c2;
      if (r === 255 && g === 255 && b === 255) { out[o] = out[o + 1] = out[o + 2] = out[o + 3] = 255; continue; }  // beyond the grid
      const cls = g >> 6, z = ((g & 63) << 8 | b) - R.zoffset, hs = r / 254 + R.hsLo;
      let ramp, t, strength = 0.9, white = 0;
      if (mode === 'bed') {
        if (z < 0) { ramp = T.bedLo; t = z / bed[0]; } else { ramp = T.bedHi; t = z / bed[1]; }
        if (cls === 0) white = 0.25;
      } else if (cls === 0) { ramp = T.ocean; t = (z - ocean[0]) / (ocean[1] - ocean[0]); strength = 0.25; }
      else if (cls === 1) { ramp = T.rock; t = z / 2000; strength = 1; }
      else if (cls === 2) { ramp = T.ice; t = (z - srf[0]) / (srf[1] - srf[0]); }
      else { ramp = null; strength = 0.6; }
      if (ramp) { const k = 3 * Math.round(n * clamp01(t)); c0 = ramp[k]; c1 = ramp[k + 1]; c2 = ramp[k + 2]; }
      else { c0 = shelf[0]; c1 = shelf[1]; c2 = shelf[2]; }
      const f = 1 + strength * hs;
      c0 = Math.min(c0 * f, 255); c1 = Math.min(c1 * f, 255); c2 = Math.min(c2 * f, 255);
      out[o] = c0 + white * (255 - c0); out[o + 1] = c1 + white * (255 - c1); out[o + 2] = c2 + white * (255 - c2);
      out[o + 3] = 255;
    }
    return out;
  }

  // 1, 2 or 5 times a power of ten giving about five steps over span
  function niceStep(span) {
    const raw = Math.max(span, 1) / 5, mag = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / mag;
    return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
  }
  function niceTicks(lo, hi) {
    const step = niceStep(hi - lo), out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(Math.round(v));
    return out;
  }
  function percentile(h, total, q) {
    let acc = 0; const goal = q * total;
    for (let i = 0; i < h.length; i++) { acc += h[i]; if (acc >= goal) return i; }
    return h.length - 1;
  }

  function init(cfg) {
    proj4.defs(cfg.epsg, cfg.proj4);
    ol.proj.proj4.register(proj4);
    const proj = ol.proj.get(cfg.epsg);
    proj.setExtent(cfg.extent);
    const grid = new ol.tilegrid.TileGrid({extent: cfg.extent, origin: cfg.origin,
                                           resolutions: cfg.resolutions, tileSize: cfg.tileSize});

    // adaptive styles: packed data tiles coloured in the browser with the limits of the view
    const R = Object.assign({shelfRgb: hexRgb(cfg.raster.shelf)}, cfg.raster), T = rampTables(R.ramps);
    const defaults = R.limits;
    const limits = cfg.styles.map(() => Object.assign({}, defaults));      // current limits per style
    const rawTiles = new Map();                                            // decoded packed tiles
    const paintBuf = new Uint8ClampedArray(4 * TILE_PX * TILE_PX);
    const adaptiveOn = () => adaptToggle.checked;

    function rawTile(i, z, x, y) {
      const key = `${i}/${z}/${x}/${y}`;
      if (!rawTiles.has(key)) {
        rawTiles.set(key, fetch(`${base}tiles/${cfg.styles[i].id}/${z}/${x}/${y}.png`)
          .then(r => { if (!r.ok) throw new Error(r.status); return r.blob(); })
          .then(b => createImageBitmap(b, {colorSpaceConversion: 'none', premultiplyAlpha: 'none'}))
          .then(bmp => {
            const ctx = new OffscreenCanvas(TILE_PX, TILE_PX).getContext('2d', {willReadFrequently: true});
            ctx.drawImage(bmp, 0, 0);
            return ctx.getImageData(0, 0, TILE_PX, TILE_PX).data;
          }).catch(e => { rawTiles.delete(key); throw e; }));
        if (rawTiles.size > RAW_CACHE) rawTiles.delete(rawTiles.keys().next().value);
      }
      return rawTiles.get(key);
    }

    function adaptiveSource(i, mode) {
      return new ol.source.DataTile({projection: proj, tileGrid: grid, transition: 0, interpolate: true,
        loader: (z, x, y) => rawTile(i, z, x, y).then(raw => {
          const cv = document.createElement('canvas'); cv.width = cv.height = TILE_PX;
          const px = paintTile(raw, mode, limits[i], R, T, paintBuf);
          cv.getContext('2d').putImageData(new ImageData(px, TILE_PX, TILE_PX), 0, 0);
          return cv;
        })});
    }

    const baseLayers = cfg.styles.map((s, i) => new ol.layer.Tile({
      visible: i === 0,
      source: s.adaptive ? adaptiveSource(i, s.adaptive)
        : new ol.source.XYZ({projection: proj, tileGrid: grid, transition: 0,
                             url: base + 'tiles/' + s.id + '/{z}/{x}/{y}.jpg'}),
    }));

    const gj = new ol.format.GeoJSON({dataProjection: proj, featureProjection: proj});
    const vector = (file, style, opts) => new ol.layer.Vector(Object.assign({
      source: new ol.source.Vector({url: base + file, format: gj}), style}, opts || {}));

    const overlays = {
      contours: vector('contours.geojson', lineStyle('rgba(40,40,40,0.45)', 0.8)),
      margin: vector('margin.geojson', f => f.get('kind') === 'grounding line'
        ? lineStyle('rgba(60,60,60,0.8)', 1) : lineStyle('rgba(0,0,0,0.75)', 1.3)),
      divides: vector('divides.geojson', lineStyle('rgba(0,0,0,0.6)', 2)),
      seaice: new ol.layer.Group({layers: cfg.seaice.map(s =>
        vector(s.file, lineStyle('rgba(29,78,137,0.8)', 1.6, s.dash)))}),
      regions: vector('regions.geojson', f => new ol.style.Style({text: new ol.style.Text({
        text: f.get('name'), font: 'bold 22px sans-serif', fill: new ol.style.Fill({color: 'rgba(40,40,40,0.8)'}),
        stroke: new ol.style.Stroke({color: 'rgba(255,255,255,0.8)', width: 4})})}), {declutter: 'regions'}),
      names: vector('labels.geojson', (f, res) => f.get('tier') > 1 && res > TIER2_RES ? null : labelStyle(f),
                    {declutter: 'names', names: true}),
      gazetteer: vector('gazetteer.geojson', labelStyle, {declutter: 'names', names: true, maxResolution: GAZ_RES}),
      graticule: new ol.layer.Graticule({visible: false, wrapX: false, showLabels: true,
        strokeStyle: new ol.style.Stroke({color: 'rgba(0,0,0,0.3)', width: 1, lineDash: [4, 4]})}),
    };

    const view = new ol.View({projection: proj, extent: cfg.extent, constrainOnlyCenter: true,
                              maxResolution: cfg.resolutions[0] * 1.5,
                              minResolution: cfg.resolutions[cfg.resolutions.length - 1] / 4});
    const popup = document.createElement('div'); popup.className = 'cryo-popup';
    const overlay = new ol.Overlay({element: popup, positioning: 'bottom-center', offset: [0, -10], autoPan: true});
    const map = new ol.Map({
      target: el, layers: [...baseLayers, ...Object.values(overlays)], view, overlays: [overlay],
      controls: ol.control.defaults.defaults({attribution: false}).extend([
        new ol.control.ScaleLine({units: 'metric'}),
        new ol.control.MousePosition({projection: 'EPSG:4326', className: 'cryo-mouse',
          coordinateFormat: c => `${fmt(Math.abs(c[1]), 2)}°${c[1] >= 0 ? 'N' : 'S'}, ` +
                                 `${fmt(Math.abs(c[0]), 2)}°${c[0] >= 0 ? 'E' : 'W'}`}),
      ]),
    });
    // keep clear of the control panel on wide screens
    const fitPadding = () => [10, el.clientWidth > 900 ? 320 : 10, 10, 10];
    view.fit(cfg.extent, {padding: fitPadding()});
    el.olMap = map;                      // for debugging from the console

    // --- popup with name details
    function showPopup(f) {
      const p = f.getProperties();
      popup.innerHTML = `<b>${p.name}</b><br><span>${p.type}</span>` +
        (p.alt ? `<br><i>also: ${p.alt}</i>` : '') +
        `<br><small>${fmt(Math.abs(p.lat), 2)}°${p.lat >= 0 ? 'N' : 'S'}, ${fmt(Math.abs(p.lon), 2)}°${p.lon >= 0 ? 'E' : 'W'}` +
        ` · ${p.source}</small>`;
      overlay.setPosition(f.getGeometry().getCoordinates());
    }
    map.on('singleclick', e => {
      const f = map.forEachFeatureAtPixel(e.pixel, f => f, {layerFilter: l => l.get('names'), hitTolerance: 4});
      f ? showPopup(f) : overlay.setPosition(undefined);
    });

    // --- control panel
    const panel = document.createElement('div'); panel.className = 'cryo-panel';
    el.appendChild(panel);
    const n = cfg.numbers;
    panel.innerHTML = `
      <input class="cryo-search" list="cryo-names" placeholder="Search a name …">
      <datalist id="cryo-names"></datalist>
      <details open><summary>Map</summary>
        ${cfg.styles.map((s, i) => `<label><input type="radio" name="style" value="${i}" ${i === 0 ? 'checked' : ''}> ${s.name}</label>`).join('')}
      </details>
      <details open><summary>Layers</summary>
        ${[['names', 'Place names', true], ['gazetteer', 'All glacier names (zoom in)', true],
           ['regions', 'Drainage regions', true], ['divides', 'Drainage divides', true],
           ['contours', '500 m surface contours', true], ['margin', 'Ice margin, grounding line', true],
           ['seaice', 'Sea-ice edge', true], ['graticule', 'Graticule', false]].map(([k, t, on]) =>
          `<label><input type="checkbox" data-layer="${k}" ${on ? 'checked' : ''}> ${t}</label>`).join('')}
      </details>
      <details open><summary>Legend</summary>
        <label class="cryo-adapt" title="Stretch the colours over the elevations in view"><input type="checkbox" checked> Adapt colours to view</label>
        <div class="cryo-legend"></div></details>
      <details><summary>Key numbers</summary>${numbersTable(cfg.numbers)}</details>`;

    // on small screens start with the panel folded, except the search box
    if (el.clientWidth < 700) panel.querySelectorAll('details').forEach(d => d.open = false);

    const adaptToggle = panel.querySelector('.cryo-adapt input');
    let active = 0;                                  // index of the visible style

    // colour bars of an adaptive style for its current limits (after `colorbars` in scripts/web.jl)
    function adaptiveBars(i) {
      const L = limits[i], sample = (ramp, f, n = 16) =>
        Array.from({length: n}, (_, k) => { const j = 3 * Math.round(255 * clamp01(f(k / (n - 1)))); return `rgb(${ramp[j]},${ramp[j + 1]},${ramp[j + 2]})`; });
      const pos = (v, [lo, hi]) => (v - lo) / (hi - lo);
      if (cfg.styles[i].adaptive === 'bed') {
        const [lo, hi] = L.bed, stops = Array.from({length: 24}, (_, k) => {
          const v = lo + (hi - lo) * k / 23, ramp = v < 0 ? T.bedLo : T.bedHi, j = 3 * Math.round(255 * clamp01(v < 0 ? v / lo : v / hi));
          return `rgb(${ramp[j]},${ramp[j + 1]},${ramp[j + 2]})`; });
        return [{label: 'Bed elevation [m]', stops, ticks: niceTicks(lo, hi).map(v => [pos(v, L.bed), String(v)])}];
      }
      const [olo, ohi] = L.ocean;
      return [{label: 'Ice-surface elevation [m]', stops: sample(T.ice, t => t),
               ticks: niceTicks(...L.srf).map(v => [pos(v, L.srf), String(v)])},
              // ocean depth runs from the shallow end (right colour of the map ramp) to the deep end
              {label: 'Ocean depth [m]', stops: sample(T.ocean, t => 1 - t),
               ticks: niceTicks(-ohi, -olo).map(d => [(d + ohi) / (ohi - olo), String(d)])}];
    }

    function legend(i) {
      const adapt = !!cfg.styles[i].adaptive && adaptiveOn();
      adaptToggle.parentNode.style.display = cfg.styles[i].adaptive ? '' : 'none';
      const bars = (adapt ? adaptiveBars(i) : cfg.styles[i].colorbars).map(b => `
        <div class="cryo-cb-label">${b.label}</div>
        <div class="cryo-cb" style="background: linear-gradient(to right, ${b.stops.join(',')})"></div>
        <div class="cryo-cb-ticks">${b.ticks.map(([p, t]) => `<span style="left:${100 * p}%">${t}</span>`).join('')}</div>`).join('');
      const ice = cfg.seaice.map(s => `<div class="cryo-line"><svg width="36" height="8"><line x1="0" y1="4" x2="36" y2="4"
        stroke="rgba(29,78,137,0.8)" stroke-width="1.6" stroke-dasharray="${s.dash.join(' ')}"/></svg> ${s.label}</div>`).join('');
      panel.querySelector('.cryo-legend').innerHTML = bars + ice;
    }

    // limits of style i for the part of the map in view: 1st to 99th percentile of the elevations
    // (surface style: of the grounded ice and, separately, the ocean), rounded outwards to nice steps
    async function viewLimits(i) {
      const mode = cfg.styles[i].adaptive, res = view.getResolution(), z = grid.getZForResolution(res);
      const size = map.getSize(), w = size[0] - (el.clientWidth > 900 ? 320 : 0);
      const a = map.getCoordinateFromPixel([0, 0]), b = map.getCoordinateFromPixel([Math.max(w, 50), size[1]]);
      if (!a || !b) return null;
      const ext = [a[0], b[1], b[0], a[1]];
      const rng = grid.getTileRangeForExtentAndZ(ext, z), full = grid.getFullTileRange(z), tres = grid.getResolution(z);
      const hists = [new Int32Array(HIST_N), new Int32Array(HIST_N)], counts = [0, 0];
      const jobs = [];
      for (let x = Math.max(rng.minX, full.minX); x <= Math.min(rng.maxX, full.maxX); x++)
        for (let y = Math.max(rng.minY, full.minY); y <= Math.min(rng.maxY, full.maxY); y++)
        jobs.push(rawTile(i, z, x, y).then(raw => {
          const te = grid.getTileCoordExtent([z, x, y]);
          const c0 = Math.max(0, Math.floor((ext[0] - te[0]) / tres)), c1 = Math.min(TILE_PX, Math.ceil((ext[2] - te[0]) / tres));
          const r0 = Math.max(0, Math.floor((te[3] - ext[3]) / tres)), r1 = Math.min(TILE_PX, Math.ceil((te[3] - ext[1]) / tres));
          for (let r = r0; r < r1; r++) for (let c = c0; c < c1; c++) {
            const o = 4 * (r * TILE_PX + c), g = raw[o + 1];
            if (raw[o] === 255 && g === 255 && raw[o + 2] === 255) continue;
            const cls = g >> 6, k = mode === 'bed' ? 0 : cls === 2 ? 0 : cls === 0 ? 1 : -1;
            if (k < 0) continue;
            hists[k][(g & 63) << 8 | raw[o + 2]]++; counts[k]++;
          }
        }));
      await Promise.all(jobs);
      const span = k => {
        const lo = percentile(hists[k], counts[k], 0.01) - R.zoffset, hi = percentile(hists[k], counts[k], 0.99) - R.zoffset;
        const step = niceStep(hi - lo);
        return [Math.floor(lo / step) * step, Math.ceil(hi / step) * step];
      };
      const out = Object.assign({}, defaults), MIN = 50;       // too few pixels: keep the scale of the whole map
      if (mode === 'bed') {
        if (counts[0] >= MIN) { const [lo, hi] = span(0); out.bed = [Math.min(lo, -50), Math.max(hi, 50)]; }
      } else {
        if (counts[0] >= MIN) { const [lo, hi] = span(0); out.srf = [lo, hi > lo ? hi : lo + 100]; }
        if (counts[1] >= MIN) { const [lo, hi] = span(1); out.ocean = [lo, Math.min(hi, 0) > lo ? Math.min(hi, 0) : lo + 100]; }
      }
      return out;
    }

    // keep the limits of the active style in step with the view (throttled while the user pans and zooms)
    let timer = null, busy = false, again = false;
    async function adapt() {
      const i = active;
      if (!cfg.styles[i].adaptive) return;
      if (busy) { again = true; return; }
      busy = true;
      try {
        const lim = adaptiveOn() ? await viewLimits(i) : Object.assign({}, defaults);
        if (lim && i === active && JSON.stringify(lim) !== JSON.stringify(limits[i])) {
          limits[i] = lim; baseLayers[i].getSource().refresh(); legend(i);
        }
      } catch (e) { console.warn('colour scale', e); }
      busy = false;
      if (again) { again = false; schedule(); }
    }
    function schedule() { if (!timer) timer = setTimeout(() => { timer = null; adapt(); }, ADAPT_MS); }
    view.on(['change:center', 'change:resolution'], schedule);
    map.once('rendercomplete', schedule);
    adaptToggle.addEventListener('change', () => { legend(active); adapt(); });

    legend(0);
    panel.querySelectorAll('input[name=style]').forEach(r => r.addEventListener('change', () => {
      active = +r.value; baseLayers.forEach((l, k) => l.setVisible(k === active)); legend(active); schedule();
    }));
    panel.querySelectorAll('input[data-layer]').forEach(c => c.addEventListener('change', () =>
      overlays[c.dataset.layer].setVisible(c.checked)));

    // --- name search over the label CSV, the gazetteer and the drainage regions
    const index = [];
    Promise.all(['labels.geojson', 'gazetteer.geojson', 'regions.geojson'].map(f => fetch(base + f).then(r => r.json()))).then(cs => {
      const names = new Map();           // name => option label (regions only)
      cs.forEach(c => c.features.forEach(ft => {
        const p = ft.properties;
        [p.name, ...(p.alt ? p.alt.split(', ') : [])].forEach(nm => index.push({key: nm.toLowerCase(), ft}));
        names.set(p.name, p.type === 'drainage region' ? `${p.alt || p.name}, drainage region` : '');
      }));
      panel.querySelector('#cryo-names').innerHTML = [...names].sort(([a], [b]) => a.localeCompare(b))
        .map(([s, l]) => `<option value="${s}"${l ? ` label="${l}"` : ''}>`).join('');
    });
    panel.querySelector('.cryo-search').addEventListener('change', e => {
      const q = e.target.value.trim().toLowerCase();
      const hit = index.find(h => h.key === q) || index.find(h => h.key.includes(q));
      if (!hit) return;
      overlay.setPosition(undefined);
      // regions: zoom to their extent
      if (hit.ft.properties.extent) return view.fit(hit.ft.properties.extent, {padding: fitPadding(), duration: 800});
      const f = gj.readFeature(hit.ft);
      view.animate({center: f.getGeometry().getCoordinates(), resolution: Math.min(view.getResolution(), 400), duration: 800},
                   () => showPopup(f));
    });
  }
})();
