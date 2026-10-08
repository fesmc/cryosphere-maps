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

  function init(cfg) {
    proj4.defs(cfg.epsg, cfg.proj4);
    ol.proj.proj4.register(proj4);
    const proj = ol.proj.get(cfg.epsg);
    proj.setExtent(cfg.extent);
    const grid = new ol.tilegrid.TileGrid({extent: cfg.extent, origin: cfg.origin,
                                           resolutions: cfg.resolutions, tileSize: cfg.tileSize});

    const baseLayers = cfg.styles.map((s, i) => new ol.layer.Tile({
      visible: i === 0,
      source: new ol.source.XYZ({projection: proj, tileGrid: grid, transition: 0,
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
        stroke: new ol.style.Stroke({color: 'rgba(255,255,255,0.8)', width: 4})})}), {declutter: 'names'}),
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
    view.fit(cfg.extent, {padding: [10, el.clientWidth > 900 ? 320 : 10, 10, 10]});
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
      <details open><summary>Legend</summary><div class="cryo-legend"></div></details>
      <details><summary>Key numbers</summary>${numbersTable(cfg.numbers)}</details>`;

    // on small screens start with the panel folded, except the search box
    if (el.clientWidth < 700) panel.querySelectorAll('details').forEach(d => d.open = false);

    function legend(i) {
      const bars = cfg.styles[i].colorbars.map(b => `
        <div class="cryo-cb-label">${b.label}</div>
        <div class="cryo-cb" style="background: linear-gradient(to right, ${b.stops.join(',')})"></div>
        <div class="cryo-cb-ticks">${b.ticks.map(([p, t]) => `<span style="left:${100 * p}%">${t}</span>`).join('')}</div>`).join('');
      const ice = cfg.seaice.map(s => `<div class="cryo-line"><svg width="36" height="8"><line x1="0" y1="4" x2="36" y2="4"
        stroke="rgba(29,78,137,0.8)" stroke-width="1.6" stroke-dasharray="${s.dash.join(' ')}"/></svg> ${s.label}</div>`).join('');
      panel.querySelector('.cryo-legend').innerHTML = bars + ice;
    }
    legend(0);
    panel.querySelectorAll('input[name=style]').forEach(r => r.addEventListener('change', () => {
      const i = +r.value; baseLayers.forEach((l, k) => l.setVisible(k === i)); legend(i);
    }));
    panel.querySelectorAll('input[data-layer]').forEach(c => c.addEventListener('change', () =>
      overlays[c.dataset.layer].setVisible(c.checked)));

    // --- name search over the label CSV and the gazetteer
    const index = [];
    Promise.all(['labels.geojson', 'gazetteer.geojson'].map(f => fetch(base + f).then(r => r.json()))).then(cs => {
      const names = new Set();
      cs.forEach(c => c.features.forEach(ft => {
        const p = ft.properties;
        [p.name, ...(p.alt ? p.alt.split(', ') : [])].forEach(nm => index.push({key: nm.toLowerCase(), ft}));
        names.add(p.name);
      }));
      panel.querySelector('#cryo-names').innerHTML = [...names].sort().map(s => `<option value="${s}">`).join('');
    });
    panel.querySelector('.cryo-search').addEventListener('change', e => {
      const q = e.target.value.trim().toLowerCase();
      const hit = index.find(h => h.key === q) || index.find(h => h.key.includes(q));
      if (!hit) return;
      const f = gj.readFeature(hit.ft);
      overlay.setPosition(undefined);
      view.animate({center: f.getGeometry().getCoordinates(), resolution: Math.min(view.getResolution(), 400), duration: 800},
                   () => showPopup(f));
    });
  }
})();
