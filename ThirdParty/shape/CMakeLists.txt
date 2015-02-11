cmake_minimum_required(VERSION 2.8.8)

project(Shapelib)

add_library(shp
	shpopen.c dbfopen.c safileio.c shptree.c
)
