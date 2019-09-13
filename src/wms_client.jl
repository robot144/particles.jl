
using HTTP
using Base64
using Dates
using Plots
using FileIO
using ImageMagick
#using TestImages

"""
Info about some WMS servers. Some of the info could be collected using the GetCapabilities request
to the server. For now this was done manually.
"""
KnownServers=Dict()
#emodnet-bathymetry
#figure out more details (manually for now)
# https://ows.emodnet-bathymetry.eu/wms/service?SERVICE=WMS&VERSION=1.1.1&request=getcapabilities
emo=Dict("scheme" => "https", 
	 "host" => "ows.emodnet-bathymetry.eu",
	 "path" => "/wms/service",
	 "service" => "WMS",
	 "version" => "1.1.1",
	 "layers" => ["emodnet:mean_atlas_land","emodnet:mean_multicolour","emodnet:mean_rainbowcolour"],
	 "styles" => ["mean_atlas_land","mean_multicolur","mean_rainbowcolour"],
	 "default_layer" => 3,
	 "srs" => "epsg:4326",
	 "format" => "image/png",
	 "extent" => [-36.0,15.0,43.0,90.0])
KnownServers["emodnet-bathymetry"]=emo
#open streetmap OSM
# http://ows.terrestris.de/osm/service?SERVICE=WMS&VERSION=1.1.1&request=getcapabilities
osm=Dict("scheme" => "https", 
	 "host" => "ows.terrestris.de",
	 "path" => "/osm/service",
	 "service" => "WMS",
	 "version" => "1.1.1",
	 "layers" => ["OSM-WMS","TOPO-WMS","TOPO-OSM-WMS","SRTM30-Colored-Hillshade"],
	 "styles" => ["default","default","default","default"],
	 "default_layer" => 1,
	 "srs" => "epsg:4326",
	 "format" => "image/png",
	 "extent" => [-180.0,-90.0,180.0,90.0])
KnownServers["open-streetmap"]=osm
# GEBCO https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv?
# polar view  https://www.gebco.net/data_and_products/gebco_web_services/north_polar_view_wms/mapserv?
gebco=Dict("scheme" => "https", 
	 "host" => "www.gebco.net",
	 "path" => "/data_and_products/gebco_web_services/web_map_service/mapserv",
	 "service" => "WMS",
	 "version" => "1.1.1",
	 "layers" => ["GEBCO_LATEST"],
	 "styles" => ["default"],
	 "default_layer" => 1,
	 "srs" => "epsg:4326",
	 "format" => "image/png",
	 "extent" => [-180.0,-90.0,180.0,90.0])
KnownServers["gebco"]=gebco

mutable struct WmsServer
   scheme::String
   host::String
   path::String
   service::String
   version::String
   layer::String #selected layer
   style::String 
   format::String
   srs::String
   extent::Array{Float64,1}
   verbose::Int64
   """ create """
   function WmsServer() #provide a default map for convenience
      label="emodnet-bathymetry"
      return WmsServer(label)
   end
   """ 
   s=WmsServer("open-streetmap",2)
   Create new object representing the WmsServer locally. Labels are from the keys in the KnownServers.
   layer_number<1 denotes the default layer.
   """
   function WmsServer(label::String,layer_number=-1)
      sdata=KnownServers[label]
      scheme=sdata["scheme"]
      host=sdata["host"]
      path=sdata["path"]
      service=sdata["service"]
      version=sdata["version"]
      ilayer=sdata["default_layer"]
      if layer_number>0
         ilayer=layer_number
      end
      layer=sdata["layers"][ilayer]
      style=sdata["styles"][ilayer]
      format=sdata["format"]
      srs=sdata["srs"]
      extent=copy(sdata["extent"])
      verbose=1
      return new(scheme,host,path,service,version,layer,style,format,srs,extent,verbose)
   end
end
"""
   img=get_map(wms_server,bbox)

Download a wms map for a given region and load it as an image.
bbox contains xmin ymin xmax ymax
"""
function get_map(wms_server,boundingbox,width=1200,height=800)
   query=Dict{String,String}()
   query["SERVICE"]=wms_server.service
   query["VERSION"]=wms_server.version
   query["REQUEST"]="GetMap"
   query["width"]="$(width)"
   query["height"]="$(height)"
   query["layers"]=wms_server.layer
   query["styles"]=wms_server.style
   query["format"]=wms_server.format
   query["SRS"]=wms_server.srs
   b=boundingbox
   query["bbox"]="$(b[1]),$(b[2]),$(b[3]),$(b[4])"
   #println(query)
   url=HTTP.URI(; scheme=wms_server.scheme, host=wms_server.host, path=wms_server.path, query=query)
   println(url)
   #handle cache
   if !isdir(".cache")
      println("Creating cache folder .cache")
      mkdir(".cache")
   else
      println("Cache folder .cache already exists")
   end
   url_hash=hash(url)
   filename=joinpath(".cache","$(wms_server.host)_$(url_hash).png")
   if !isfile(filename)
      println("Cache file $(filename) does not exist yet. Download new map.")
      HTTP.download(string(url),filename ; verbose=wms_server.verbose)
   else
      println("Cache file $(filename) already exists. Using cached file")
   end
   img = load(filename)
   return img
end


"""
   plot_image(image,bbox)

Plot and image scaled to a given range, such that other elements added to the graph 
will appear in the right position.
image
bbox contains xmin ymin xmax ymax
"""
function plot_image(image,boundingbox)
   bbox=boundingbox
   f=plot([bbox[1], bbox[3]], [bbox[2], bbox[4]], image[end:-1:1, :], yflip = false)
   return f
end

