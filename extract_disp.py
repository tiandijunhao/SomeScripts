i=0
j=0
k=0
num_atom=0
position0=[]
position=[]
disp=[]
disp_f=[]
latticetmp1=[]
latticetmp2=[]
latticetmp3=[]
s_step=int(input("starting step:\n"))
d_step=int(input("interval of steps:\n"))

with open('PPOSCAR','r') as pposcar:
	text=pposcar.readlines()
	coefficient=float(text[1].strip())
	for a in range(3):
		latticetmp1.append(float(text[2].strip().split('   ')[a])*coefficient)
	for a in range(3):
		latticetmp2.append(float(text[3].strip().split('   ')[a])*coefficient)
	for a in range(3):
		latticetmp3.append(float(text[4].strip().split('   ')[a])*coefficient)
	#print(latticetmp1,latticetmp2,latticetmp3)
	if text[5].strip() < text[6].strip():
		print('please check the input format of PPOSCAR')
	else:
		num_class=len(text[5].strip().split(' '))
		for c in range(num_class):
			num_atom=num_atom+int(text[6].strip().split(' ')[c])
		for n in range(num_atom):
			position0.append(text[n+8].split())
# print(num_atom)


with open('vasprun.xml','r') as	vasprun:
	text=vasprun.readlines()
	num_lines=len(text)
	for n in range(num_lines):
		
		# # !!! old the original positions of atoms
		# if text[n] == "  <varray name=\"positions\" >\n":
			# for m in range(num_atom):
				# position0.append(text[n+m+1].split())
		
		if text[n] == "   <varray name=\"positions\" >\n":
			i=i+1
			if i > s_step-1:
				j=j+1
				if ((j-1)-d_step*((j-1)//d_step))==0:
					k=k+1
					#position.append(i)
					for m in range(num_atom):
						position.append(text[n+m+1].split())



print("total_steps=",i)
#print(j)
print("output_steps=",k)

for each_line in position0:
	for each_item in each_line:
		each_item = float(each_item)


# with open('0','w') as f0:
	# for n in range(len(position0)):
		# for m in range(2):
			# print(position0[n][m],end=' ', file=f0)
		# print(position0[n][2],end='\n', file=f0)

for each_line in position:
	each_line.pop(4)
	each_line.pop(0)
	for each_item in each_line:
		each_item = float(each_item)

# with open('1','w') as f1:
	# for n in range(len(position)):
		# for m in range(2):
			# print(position[n][m],end=' ', file=f1)
		# print(position[n][2],end='\n', file=f1)


# differences of positions equal disps
#print(len(position),num_atom)
time_step=int(len(position)/num_atom)
for n in range(time_step):
	for m in range(num_atom):
		for a in range(3):
			disp.append('%.8f' %(float(position[m+n*num_atom][a])-float(position0[m][a])))
			
for n in range(int(len(disp)/3)):
	if float(disp[n*3]) > 0.5:
		disp[n*3]=float(disp[n*3])-1
	elif float(disp[n*3]) < -0.5:
		disp[n*3]=float(disp[n*3])+1
	if float(disp[n*3+1]) > 0.5:
		disp[n*3+1]=float(disp[n*3+1])-1
	elif float(disp[n*3+1]) < -0.5:
		disp[n*3+1]=float(disp[n*3+1])+1
	if float(disp[n*3+2]) > 0.5:
		disp[n*3+2]=float(disp[n*3+2])-1
	elif float(disp[n*3+2]) < -0.5:
		disp[n*3+2]=float(disp[n*3+2])+1

for n in range(int(len(disp)/3)):
	disp_f.append(float(disp[3*n])*latticetmp1[0]+float(disp[3*n+1])*latticetmp2[0]+float(disp[3*n+2])*latticetmp3[0])
	disp_f.append(float(disp[3*n])*latticetmp1[1]+float(disp[3*n+1])*latticetmp2[1]+float(disp[3*n+2])*latticetmp3[1])
	disp_f.append(float(disp[3*n])*latticetmp1[2]+float(disp[3*n+1])*latticetmp2[2]+float(disp[3*n+2])*latticetmp3[2])

for n in range(len(disp_f)):
	disp_f[n]=float(disp_f[n])/0.5291772

with open('disp.dat','w') as fd:
	for n in range(int(len(disp)/3)):
		for m in range(2):
			if float(disp_f[n*3+m]) > 0:
				print('     ','%.8f' % float(disp_f[n*3+m]), end=' ', file=fd)
			else:
				print('    ','%.8f' % float(disp_f[n*3+m]), end=' ', file=fd)
		if float(disp_f[n*3+2]) > 0:
			print('     ','%.8f' % float(disp_f[n*3+2]), end='\n', file=fd)
		else:
			print('    ','%.8f' % float(disp_f[n*3+2]), end='\n', file=fd)


