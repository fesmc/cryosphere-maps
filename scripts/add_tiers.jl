# One-off: append a `tier` column (1 = major, 2 = detail) to the label gazetteers.
using DelimitedFiles

const TIER2 = Dict(
"antarctica" => Set(["Larsen D Ice Shelf","Nickerson Ice Shelf","Crosson Ice Shelf","Dotson Ice Shelf","Cook Ice Shelf",
  "Mertz Glacier Tongue","Ninnis Glacier Tongue","Nansen Ice Shelf","Jelbart Ice Shelf","Ekström Ice Shelf",
  "Lazarev Ice Shelf","Roi Baudouin Ice Shelf","Stancomb-Wills Glacier Tongue","Pine Island Ice Shelf","Thwaites Ice Shelf",
  "Smith Glacier","Pope Glacier","Kohler Glacier","Nimrod Glacier","Shackleton Glacier","Axel Heiberg Glacier","Scott Glacier",
  "Amundsen Glacier","Mawson Glacier","Mercer Ice Stream","Bindschadler Ice Stream","MacAyeal Ice Stream","Evans Ice Stream",
  "Carlson Inlet","Academy Glacier","Support Force Glacier","Bailey Ice Stream","Stancomb-Wills Glacier","Mellor Glacier",
  "Fisher Glacier","Frost Glacier","Holmes Glacier","Dibble Glacier","Mertz Glacier","Ninnis Glacier","Cook Glacier",
  "Rennick Glacier","Lillie Glacier","Hays Glacier","Rayner Glacier","Kemp Land","Oates Land","Kunlun Station",
  "Dome Fuji Station","Concordia Station","Neumayer III","Halley VI","Casey","Davis","Mawson","Dumont d'Urville","Rothera",
  "McMurdo Station","Law Dome","Talos Dome","Siple Dome","Davis Sea","Cosmonaut Sea","Riiser-Larsen Sea","Lazarev Sea",
  "Cooperation Sea","Dumont d'Urville Sea","Somov Sea","Mawson Sea","Scotia Sea","Siple Island","Adelaide Island",
  "Ross Island","Thurston Island","Ellsworth Land","Antarctic Peninsula"]),
"greenland" => Set(["C.H. Ostenfeld Gletscher","Tracy Gletscher","Hayes Gletscher","Steenstrup Gletscher","Køge Bugt Gletscher",
  "Ikertivaq Gletscher","Midgård Gletscher","Sermeq Avannarleq","Store Gletscher","Kong Oscar Gletscher","Rink Isbræ",
  "Kangiata Nunaata Sermia","Renland","Summit Station","Danmarkshavn","Station Nord","Upernavik","Narsarsuaq",
  "Ittoqqortoormiit","Nares Strait","Lincoln Sea","Disko Bay","Melville Bay","Germania Land","Inglefield Land",
  "Washington Land","Jameson Land","Kangerlussuaq Fjord","Ilulissat Icefjord","Nuup Kangerlua (Godthåbsfjord)"]),
)

for reg in keys(TIER2)
    f = joinpath(@__DIR__, "..", "data", "labels_$(reg).csv")
    d, h = readdlm(f, ',', Any; header=true, quotes=true)
    "tier" in h && continue
    names = strip.(string.(d[:, 1]))
    missing_ = setdiff(TIER2[reg], names); isempty(missing_) || @warn "not found" reg missing_
    tier = [n in TIER2[reg] ? 2 : 1 for n in names]
    q(s) = occursin(',', string(s)) ? "\"$(s)\"" : string(s)
    open(f, "w") do io
        println(io, join(vcat(vec(h), "tier"), ','))
        for (k, r) in enumerate(eachrow(d))
            println(io, join(vcat(q.(collect(r)), tier[k]), ','))
        end
    end
    println(reg, ": tier1 = ", count(==(1), tier), ", tier2 = ", count(==(2), tier))
end
