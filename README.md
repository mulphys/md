# ReMoDy: Reactive Molecular Dynamics Solver

## Author: [Andrei V Smirnov](mailto:andrei.v.smirnov@gmail.com)

## [URL](http://galacticbubble.com/remody)

## Molecular dynamics simulation for kinetic reactions and interfacial chemistry

### SUB-DIRECTORIES:

mp          - multi-processor version with local time-stepping scheme
sp          - single processor version with standard time-stepping scheme
doc         - documents directory
util        - some utilities for pre- and post-processing

### IMPORTANT FILES:

remody.tgz                - TAR-gzipped archive of multi-processor version with local time-stepping scheme
doc/usersguide.pdf        - User's Guide (PDF)
doc/programmersguide.pdf  - Programmer's Guide (PDF)
doc/html/index.html       - Source-code documentaiton (HTML)

### RETRIEVING FROM ARCHIVE:

tar xvzf remody.tbz

### COMPILING:

- Single processor version:
	cd sp; make

- Multi processor version:
	cd mp; make


### CONFIGURING INPUT:

Edit the input XML file (see documentation in /doc).


### RUNNING:

cd run
./job

or for GUI version:

cd run
./gui -f job.xml

or:

./gui -f test.xml

depending on which case to run. 


## DOCUMENTATION:

### [HTML](http://doc/html/index.html)

### [PDF](doc)

### [User Guide](doc/usersguide.pdf)

### [Programmer's Guide](doc/programmersguide.pdf)


