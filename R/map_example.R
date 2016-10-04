library(sp)
library(rgeos)

load("~/Desktop/spydf_states.RData")

# many geos functions require projections and you're probably going to end
# up plotting this eventually so we convert it to albers before cleaning up
# the polygons since you should use that if you are plotting the US
spydf_states <- spTransform(spydf_states, 
                            CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96"))

# simplify the polgons a tad (tweak 0.00001 to your liking)
spydf_states <- gSimplify(spydf_states, tol = 0.00001)

# this is a well known R / GEOS hack (usually combined with the above) to 
# deal with "bad" polygons
spydf_states <- gBuffer(spydf_states, byid=TRUE, width=0)

# any bad polys?
sum(gIsValid(spydf_states, byid=TRUE)==FALSE)

usa <- gUnaryUnion(map2SpatialPolygons(spydf_states, IDs=spydf_states$names,
                                 proj4string=CRS("+proj=longlat +datum=WGS84")))

#plot(spydf_states)
