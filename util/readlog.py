import sys
import re
if len(sys.argv) != 2:
	sys.exit("Must provide the specie name, like:\n\t"+sys.argv[0]+" H2O\n\n \
		Or variable name, like:\n\t"+sys.argv[0]+" T=\n")
name=sys.argv[1]
namelen=len(name)
oldtime=''
var=''
num=re.compile('[-e\.\d]')
while True:
	line=sys.stdin.readline()
	if not line:
		break
	i=line.find('Time')
	if i>=0:
		i=line.find('=')+1
		j=line.find(' ',i)
		time=line[i:j]+' '
		if time == oldtime: continue
		oldtime=time
		if name.find('=')>0: # retrieve inline variable
			i=line.find(name)
			i=line.find('=',i)+1
			j=i
			while num.match(line[j:j+1]): j=j+1
			sys.stdout.write("%s%s\n"%(time,line[i:j]))
			continue;
		while True:
			line=sys.stdin.readline() # header
			i=0
			while True:
				i=line.find(name,i)
				if i>0:
					if not line[i+namelen].isspace():
						i=i+namelen
						continue
					break
				break
			if i>0: break
			if not sys.stdin.readline(): # skip values
				sys.exit('Unexpected end of file')
			# end while
		if not line:
			print 'Unexpected end of file'
			sys.exit(1)
		# Count columns:
		col=-1
		j=0
		while j<i:
			if line[j] == '%': col=col+1
			j=j+1
		if col<0:
			sys.exit('Column format error')
		line=sys.stdin.readline() # read values
		i=0
		for k in range(col):
			while not line[i].isspace(): 
				i=i+1
			while line[i].isspace():
				i=i+1
		if line[i].isdigit():
			j=i
			while not line[j].isspace():
				j=j+1
			conc=line[i:j]
			sys.stdout.write(time+"\t"+conc+"\n")
