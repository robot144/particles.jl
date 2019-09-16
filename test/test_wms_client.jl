
#using HTTP
#using Base64
#using Dates
using Plots
#using FileIO

#
# testing
#
function test1()
   #width=1200
   #height=1000
   #Plots.default(:size,[width,height])
   emodnet_server=WmsServer("emodnet-bathymetry")
   @test emodnet_server.host=="ows.emodnet-bathymetry.eu"
   @test emodnet_server.format=="image/png"
   @test emodnet_server.extent==[-36.0, 15.0, 43.0, 90.0]
   @test emodnet_server.srs=="epsg:4326"
   bbox=[0.0,49.0,9.0,55.0]
   img=get_map(emodnet_server,bbox)
   @test size(img)==(800, 1200)
   @test Float64(img[1,1].r)==0.7137254901960784 #carefull, since this is sensitive to server changes
   #@test isfile(joinpath(".cache","ows.emodnet-bathymetry.eu_8803805663758041344.png")) 
   plot_image(img,bbox)
 end
   
function test2()
   width=1200
   height=1000
   #Plots.default(:size,[width,height])
   osm_server=WmsServer("open-streetmap")
   @test osm_server.host=="ows.terrestris.de"
   bbox=[0.0,49.0,9.0,55.0]
   img=get_map(osm_server,bbox,width,height)
   @test size(img)==(1000, 1200)
   #@test isfile(joinpath(".cache","ows.terrestris.de_9092741814141544534.png")) 
   plot_image(img,bbox)
 end
   
function test3()
   width=2000
   height=1000
   #Plots.default(:size,[width,height])
   gebco_server=WmsServer("gebco")
   @test gebco_server.host=="www.gebco.net"
   bbox=[-180.0,-90.0,180.0,90.0]
   img=get_map(gebco_server,bbox,width,height)
   @test size(img)==(1000, 2000)
   #@test isfile(joinpath(".cache","www.gebco.net_7532605064285132715.png")) 
   plot_image(img,bbox)
 end

test1()
test2()
test3()
