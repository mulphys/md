import sys
if len(sys.argv)<2:
	print "Usage: "+sys.argv[0]+" filename\n"
	sys.exit(1)
filename=sys.argv[1]
isep=filename.rfind('/')+1
idot=filename.rfind('.')
if idot<0: idot=len(filename)
name=filename[isep:idot]
#name=filename.split('.')[0]+"-V2"
file=open(filename,'r')
vel2=0.0 # velocity squared
n=0
for line in file:
	vels=line.split()
	v2=0.0
	for i in range(len(vels)):
		v2+=float(vels[i])**2
	vel2+=v2
	n+=1

print "%s\t%g\n"%(name,vel2/n)
