
The program should be build into the ```build/``` directory, else the relative paths don't point to the correct locations.

Necessary for the Voronoi tesselation, download and install the following library:
- [http://math.lbl.gov/voro++/download/](http://math.lbl.gov/voro++/download/)

If you don't want gifs just ignore the following paragraph.
If gifs of the images should be created, ```ImageMagick``` needs to be installed.
On Mac ```brew install ImageMagick```, it is then called with ```convert```.
Also ```createAgif(fname);``` in ```image.c line 295``` needs to be commented in.

