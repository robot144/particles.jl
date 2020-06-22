# Some routines to handle AIS ship data
# Each track is assumed to be one csv file.
# The format should be like:
#MMSI,229068000,,,,,,,,,,Length,,,Width,,,,,,
#Date and UTC Time,Lat,Lon,SOG,HDG,COG,ROT,Status,IMO Number,Name,Call Sign,Bow,Stern,Overall,Port,Starboard,Overall,Draught,Destination,Vessel Type,Extra Info
#12/02/2019 18:33:21,49.776873,-2.986952,9.6,63,63.6,0,Under way using engine,9481702,ABML EVA,9HA3041,215,39,254,27,16,43,,,Cargo,N/A
#12/02/2019 18:33:30,49.777072,-2.986328,9.7,63,64,0,Under way using engine,9481702,ABML EVA,9HA3041,215,39,254,27,16,43,,,Cargo,N/A
#12/02/2019 18:33:41,49.777262,-2.985705,9.7,64,64.5,0,Under way using engine,9481702,ABML EVA,9HA3041,215,39,254,27,16,43,,,Cargo,N/A

using DataFrames
using Dates
using CSV

#for later data, the format is a bit different. I would suggest to make a separate reading routine
#until the format is somewhat settled
#for new data
#t=[]
#I use only the new Eastbound data here/Change absolute path
#push!(t,CSV.read("English Channel Set Oct 2019 Eastbound.txt",header=1));
#Change format of data headers
#for s=1:1
#rename!(t[s], Symbol("Latitude (°)")=>Symbol("Lat")) ;
#rename!(t[s], Symbol("Longitude (°)")=>Symbol("Lon"));
#rename!(t[s], Symbol("SOG (kts)")=>Symbol("SOG"));
#rename!(t[s], Symbol("COG (°)")=>Symbol("COG"));
#rename!(t[s], Symbol("HDG (°)")=>Symbol("HDG"));
#rename!(t[s], Symbol("ROT (°/min)")=>Symbol("ROT"));
#Choose random ship from data
#    t[1]=t[1][t[1][:MMSI] .==t[1].MMSI[18], : ];
#end


function read_ais_data(filename)
   println("Reading file: $(filename)")
   data=CSV.read(filename,header=2)
   #convert ASCII time-strings to DateTime format
   times=data[:,1]
   ts=[DateTime(t,"dd/mm/yyyy HH:MM:SS") for t=times]
   data.Times=ts
   return data
end

