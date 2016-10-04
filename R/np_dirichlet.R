library(maps)
library(maptools)
library(spatstat)
library(rgeos)

#lower48 <- map("usa", fill=TRUE, plot=FALSE)
#
#usa <- gUnaryUnion(map2SpatialPolygons(lower48, IDs=lower48$names,
#                                 proj4string=CRS("+proj=longlat +datum=WGS84")))

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

W <- as(spydf_states, "owin")

X <- ppp(x = c(-68.21,-109.57,-102.50,-103.25,-80.08,-107.72,-112.18,-109.93,-111.17,-104.44,-119.42,-80.78,-122.1,-81.55,-116.82,-82.87,-80.93,-114.00,-112.14,-110.80,-114.30,-105.51,-83.53,-104.87,-93.05,-88.55,-115.90,-118.55,-121.51,-86.10,-108.49,-121.75,-121.20,-123.50,-109.78,-121.16,-124.00,-105.58,-110.50,-118.68,-78.35,-103.45,-92.88,-103.48,-110.50,-119.50,-113.05), y = c(44.35, 38.68, 43.75, 29.25, 25.65, 38.57, 37.57, 38.20, 38.20, 32.17, 34.01, 33.78, 42.94, 41.24, 36.24, 24.63, 25.32, 48.80, 36.06, 43.73, 38.98, 37.73, 35.68, 31.92, 34.51, 48.10, 33.79, 36.80, 40.49, 37.18, 37.18, 46.85, 48.70, 47.97, 35.07, 36.48, 41.30, 40.40, 32.25, 36.43, 38.53, 46.97, 48.50, 43.57, 44.60, 37.83, 37.30), window = W)

y <- dirichlet(X) # Dirichlet tessellation
plot(y) # Plot tessellation
plot(X, add = TRUE) # Add points
tile.areas(y) #Areas

#counties <- map("county", c('minnesota,ramsey','minnesota,hennepin'),fill=TRUE, plot=FALSE)
#merged <- gUnaryUnion(map2SpatialPolygons(counties, IDs=counties$names, proj4string=CRS("+proj=longlat +datum=WGS84")))

