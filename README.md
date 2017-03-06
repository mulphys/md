# ReMoDy: Reactive Molecular Dynamics Solver

## Author: [Andrei V Smirnov](mailto:andrei.v.smirnov@gmail.com)

## [URL](http://galacticbubble.com/remody)

## Molecular dynamics simulation for kinetic reactions and interfacial chemistry

### SUB-DIRECTORIES:

| [mp](mp/)     | multi-processor version with local time-stepping scheme
| [sp](sp/)     | single processor version with standard time-stepping scheme
| [doc](doc/)   | documents directory
| [util](util/) | some utilities for pre- and post-processing

### IMPORTANT FILES:

| [remody.tgz](remody.tgz)               | TAR-gzipped archive of multi-processor version with local time-stepping scheme
| [doc/usersguide.pdf](remody.tgz)       | User's Guide (PDF)
| [doc/programmersguide.pdf](remody.tgz) | Programmer's Guide (PDF)
| [doc/html/index.html](remody.tgz)      | Source-code documentaiton (HTML)

### RETRIEVING FROM ARCHIVE:

tar xvzf remody.tbz

### COMPILING:

#### Single processor version:

cd sp; make

#### Multi processor version:

cd mp; make


### CONFIGURING INPUT:

Edit the input XML file (see documentation in /doc).


### RUNNING DEMOS:

#### GUI version:

cd run; ./gui -f syngas.xml


#### Batch job:

cd mp/run; ./job

## DOCUMENTATION:

### [HTML](http://galacticbubble.com/remody/doc/html/index.html)

### [PDF](doc)

### [User Guide](doc/usersguide.pdf)

### [Programmer's Guide](doc/programmersguide.pdf)


