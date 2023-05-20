i=0
j=0
k=0
num_atom=0
force=[]
s_step=int(input("starting step:\n"))
d_step=int(input("interval of steps:\n"))

with open('PPOSCAR','r') as pposcar:
	text=pposcar.readlines()
	if text[5].strip() < text[6].strip():
		print('please check the input format of PPOSCAR')
		# for n in range(num_atom):
			# position0.append(text[n+7].split())
	else:
		num_class=len(text[5].strip().split(' '))
		for c in range(num_class):
			num_atom=num_atom+int(text[6].strip().split(' ')[c])
		# for n in range(num_atom):
			# position0.append(text[n+8].split())
	# print(num_atom)

with open('vasprun.xml','r') as	vasprun:
	text=vasprun.readlines()
	num_lines=len(text)
	for n in range(num_lines):

		# the original positions of atoms
		# if text[n] == "  <varray name=\"positions\" >\n":
			# for m in range(num_atom):
				# position0.append(text[n+m+1].split())
		if text[n] == "  <varray name=\"forces\" >\n":
			i=i+1
			if i > s_step-1:
				j=j+1
				if ((j-1)-d_step*((j-1)//d_step))==0:
					k=k+1
					#force.append(i)
					for m in range(num_atom):
						force.append(text[n+m+1].split())

print("total_steps=",i)
#print(j)
print("output_steps=",k)

for each_line in force:
	each_line.pop(4)
	each_line.pop(0)
	for each_item in each_line:
		each_item = float(each_item)

with open('force.dat','w') as ff:
	for n in range(len(force)):
		for m in range(2):
			if float(force[n][m]) > 0:
				print('     ', '%9e' %float(float(force[n][m])/25.71103),end=' ', file=ff)
			else:
				print('    ', '%9e' %float(float(force[n][m])/25.71103),end=' ', file=ff)
		if float(force[n][2]) > 0:
			print('     ', '%9e' %float(float(force[n][2])/25.71103),end='\n', file=ff)
		else:
			print('    ', '%9e' %float(float(force[n][2])/25.71103),end='\n', file=ff)

