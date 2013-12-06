CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

PROJECT(Shapelib)

ADD_LIBRARY(shp
	shpopen.c dbfopen.c safileio.c shptree.c
)
