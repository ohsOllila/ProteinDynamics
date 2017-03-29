#!/usr/bin/python3
import mdtraj as md
import math
import numpy
import sys
#import math3d

#print(sys.argv[1], sys.argv[2])
traj = md.load(sys.argv[1], top=sys.argv[2])
inertia=md.compute_inertia_tensor(traj)
axis1=[[0]*3]*traj.n_frames
axis2=[[0]*3]*traj.n_frames
axis3=[[0]*3]*traj.n_frames
eval1=[0]*traj.n_frames
eval2=[0]*traj.n_frames
eval3=[0]*traj.n_frames
corr0=[0]*traj.n_frames
corr1=[0]*traj.n_frames
corr2=[0]*traj.n_frames
alpha1=[0]*traj.n_frames
alpha2=[0]*traj.n_frames
alpha3=[0]*traj.n_frames
sum1=[0]*traj.n_frames
lag=traj.n_frames/10
timestep=traj.time[1]-traj.time[0]

def cdot(axisI,axisIprev,axisIIprev):
  return [numpy.dot(axisIprev, axisI)/(numpy.sum(axisIprev**2)**(1/2)),numpy.dot(axisIIprev, axisI)/(numpy.sum(axisIIprev**2)**(1/2))]

#e_values=[[0]*3]*traj.n_frames
#e_vectors=[[[0]*3]*3]*traj.n_frames
EVALUESfile=open(str(sys.argv[3]) + "EVALUES.dat", 'w+')
j=0
while j<traj.n_frames:
    e_values1, e_vectors1 = numpy.linalg.eig(inertia[j])
    for i in range(len(e_values1)):
	# find biggest eigen value
        if e_values1[i]==max(e_values1):
            eval1[j]=e_values1[i]
            axis1[j]=e_vectors1[:,i]
	    # find smallest eigen value
        elif e_values1[i]==min(e_values1):
            eval3[j]=e_values1[i]
            axis3[j]=e_vectors1[:,i]
	    # middle eigen value
        else:
            eval2[j]=e_values1[i]
            axis2[j]=e_vectors1[:,i]
    #cdot1=[numpy.dot(axis1[0], axis1[j])/(numpy.sum(axis1[0]**2)**(1/2)),numpy.dot(axis2[0], axis1[j])/(numpy.sum(axis2[0]**2)**(1/2))]
    #cdot2=[numpy.dot(axis2[0], axis2[j])/(numpy.sum(axis2[0]**2)**(1/2)),numpy.dot(axis3[0], axis2[j])/(numpy.sum(axis3[0]**2)**(1/2))]
    #cdot3=[numpy.dot(axis3[0], axis3[j])/(numpy.sum(axis3[0]**2)**(1/2)),numpy.dot(axis1[0], axis3[j])/(numpy.sum(axis1[0]**2)**(1/2))]
    #alpha1[j]=math.atan(cdot1[1]/cdot1[0])
    #alpha2[j]=math.atan2(cdot2[1],cdot2[0])
    #alpha3[j]=math.atan2(cdot3[1],cdot3[0])
    if j>1:
        if numpy.dot(axis1[j-1],axis1[j])<0:
            axis1[j]=(-1)*axis1[j]
        if numpy.dot(axis2[j-1],axis2[j])<0:
            axis2[j]=(-1)*axis2[j]
        if numpy.dot(axis3[j-1],axis3[j])<0:
            axis3[j]=(-1)*axis3[j]
        #cdot1=cdot(axis1[j],axis1[j-1],axis2[j-1])
        #cdot1TST=cdot(axis2[j],axis1[j-1],axis2[j-1])
        #print(j, numpy.dot(axis1[j-1], axis1[j]),numpy.dot(axis1[j-1], axis2[j]),numpy.dot(axis3[j-1], axis3[j]))# ,(alpha1[j]-alpha1[j-1])*180/3.14, (alpha2[j]-alpha2[j-1])*180/3.14, (alpha3[j]-alpha3[j-1])*180/3.14)
        #print(j,math.sqrt(math.atan(cdot1[1]/cdot1[0])**2)*180/3.14,math.sqrt(math.atan(cdot1TST[1]/cdot1TST[0])**2)*180/3.14)
        #if math.sqrt(math.atan(cdot1[1]/cdot1[0])**2)>math.sqrt(math.atan(cdot1TST[1]/cdot1TST[0])**2):
            #tmp=axis1[j]
            #axis1[j]=axis2[j]
            #axis2[j]=tmp
            #print(j,math.atan(cdot1[1]/cdot1[0])*180/3.14,math.atan(cdot1TST[1]/cdot1TST[0])*180/3.14,numpy.dot(numpy.cross(axis1[j],axis2[j]),axis3[j]))
        #if numpy.dot(axis1[j-1],axis1[j])<0:
        #    axis1[j]=(-1)*axis1[j]
        #if numpy.dot(axis2[j-1],axis2[j])<0:
        #    axis2[j]=(-1)*axis2[j]
        #if numpy.dot(axis3[j-1],axis3[j])<0:
        #    axis3[j]=(-1)*axis3[j]
        #cdot1=cdot(axis1[j],axis1[j-1],axis2[j-1])
        #cdot1TST=cdot(axis2[j],axis1[j-1],axis2[j-1])
        #cdot2=cdot(axis2[j],axis2[j-1],axis3[j-1])
        #print(j,math.atan(cdot1[1]/cdot1[0])*180/3.14,math.atan(cdot1TST[1]/cdot1TST[0])*180/3.14,math.atan(cdot2[1]/cdot2[0])*180/3.14)

        #if numpy.dot(axis1[j-1],axis1[j])<0:
            #r = math3d.Orientation.new_axis_angle(axis2[j], math.pi)
            #v = math3d.Vector(axis1[j])
            #axis1[0]=(r * v)[0]
            #axis1[1]=(r * v)[1]
            #axis1[2]=(r * v)[2]
        #    axis1[j]=(-1)*axis1[j]
            #axis3[j]=numpy.cross(axis1[j],axis2[j])
        #if numpy.dot(axis2[j-1],axis2[j])<0:
            #r = math3d.Orientation.new_axis_angle(axis3[j], math.pi)
            #v = math3d.Vector(axis2[j])
            #axis2[0]=(r * v)[0]
            #axis2[1]=(r * v)[1]
            #axis2[2]=(r * v)[2]
        #    axis2[j]=(-1)*axis2[j]
        #if numpy.dot(axis3[j-1],axis3[j])<0:
            #r = math3d.Orientation.new_axis_angle(axis1[j], math.pi)
            #v = math3d.Vector(axis3[j])
            #axis3[0]=(r * v)[0]
            #axis3[1]=(r * v)[1]
            #axis3[2]=(r * v)[2]
        #    axis3[j]=(-1)*axis3[j]
            #    tmpEVAL=eval3[j]
            #    tmpAXIS=axis3[j]
            #    eval3[j]=eval2[j]
            #    axis3[j]=axis3[j]
            #    eval2[j]=tmpEVAL
            #    axis2[j]=tmpAXIS
            #print(numpy.cross(axis1[j],axis2[j]),axis3[j])
        #if numpy.dot(axis1[j-1],axis2[j])>numpy.dot(axis1[j-1],axis1[j]) and numpy.dot(axis1[j-1],axis1[j])<0.55:
        #    print(j, numpy.dot(axis1[j-1], axis2[j]),numpy.dot(axis1[j-1], axis1[j]))
        #    tmp=axis1[j]
        #    axis1[j]=axis2[j]
        #    axis2[j]=tmp
        #    if numpy.dot(axis1[j-1],axis1[j])<0:
        #        axis1[j]=(-1)*axis1[j]
        #    if numpy.dot(axis2[j-1],axis2[j])<0:
        #        axis2[j]=(-1)*axis2[j]
        #    print(j, numpy.dot(axis1[j-1], axis2[j]),numpy.dot(axis1[j-1], axis1[j]),numpy.dot(axis2[j-1], axis2[j]),numpy.dot(axis3[j-1], axis3[j]))
        #print(numpy.dot(axis1[j-1],axis1[j]),numpy.dot(axis2[j-1],axis2[j]),numpy.dot(axis3[j-1],axis3[j]))
        #print(j,axis1[j],axis2[j],axis3[j],numpy.dot(axis1[j-1],axis1[j]))
    print(j, eval1[j], eval2[j], eval3[j],file=EVALUESfile)
    j=j+1

ANGLESfile=open(str(sys.argv[3]) + "ANGLES.dat", 'w+')    

j=1
while j<traj.n_frames:
    #    cdot1=[numpy.dot(axis1[j-1], axis1[j])/(numpy.sum(axis1[j-1]**2)**(1/2)),numpy.dot(axis2[j-1], axis1[j])/(numpy.sum(axis2[j-1]**2)**(1/2))]
    #    cdot2=[numpy.dot(axis2[j-1], axis2[j])/(numpy.sum(axis2[j-1]**2)**(1/2)),numpy.dot(axis3[j-1], axis2[j])/(numpy.sum(axis3[j-1]**2)**(1/2))]
    #    cdot3=[numpy.dot(axis3[j-1], axis3[j])/(numpy.sum(axis3[j-1]**2)**(1/2)),numpy.dot(axis1[j-1], axis3[j])/(numpy.sum(axis1[j-1]**2)**(1/2))]
    cdot1=cdot(axis1[j],axis1[j-1],axis2[j-1])
    cdot2=cdot(axis2[j],axis2[j-1],axis3[j-1])
    cdot3=cdot(axis3[j],axis3[j-1],axis1[j-1])
    cdot1TST=cdot(axis2[j],axis1[j-1],axis2[j-1])
    cdot2TST=cdot(axis1[j],axis2[j-1],axis3[j-1])
    if math.sqrt(math.atan(cdot1[1]/cdot1[0])**2)>math.sqrt(math.atan(cdot1TST[1]/cdot1TST[0])**2):
        print("Axes swapped at time frame:",j,",  new angle:",math.atan(cdot1TST[1]/cdot1TST[0])*180/3.14)
        alpha1[j]=alpha1[j-1]+math.atan(cdot1TST[1]/cdot1TST[0])
        alpha2[j]=alpha2[j-1]+math.atan(cdot2TST[1]/cdot2TST[0])
    else:
        alpha1[j]=alpha1[j-1]+math.atan(cdot1[1]/cdot1[0])
        alpha2[j]=alpha2[j-1]+math.atan(cdot2[1]/cdot2[0])
    alpha3[j]=alpha3[j-1]+math.atan(cdot3[1]/cdot3[0])
    if math.sqrt((alpha1[j]-alpha1[j-1])**2)*180/3.14>30 or math.sqrt((alpha2[j]-alpha2[j-1])**2)*180/3.14>30 or math.sqrt((alpha3[j]-alpha3[j-1])**2)*180/3.14>30:
      print("Angle change larger than 30 degrees at frame",j,(alpha1[j]-alpha1[j-1])*180/3.14,(alpha2[j]-alpha2[j-1])*180/3.14,(alpha3[j]-alpha3[j-1])*180/3.14)
    #print(j, numpy.dot(axis1[j-1], axis1[j]),numpy.dot(axis1[j-1], axis2[j]),numpy.dot(axis3[j-1], axis3[j]))# ,(alpha1[j]-alpha1[j-1])*180/3.14, (alpha2[j]-alpha2[j-1])*180/3.14, (alpha3[j]-alpha3[j-1])*180/3.14)
    #print(j,math.sqrt(math.atan(cdot1[1]/cdot1[0])**2)*180/3.14,math.sqrt(math.atan(cdot1TST[1]/cdot1TST[0])**2)*180/3.14)
    print(j, alpha1[j]*180/3.14, alpha2[j]*180/3.14, alpha3[j]*180/3.14,file=ANGLESfile)
    #print(j, numpy.dot(axis1[j-1], axis1[j]),numpy.dot(axis2[j-1], axis2[j]),numpy.dot(axis3[j-1], axis3[j]) ,(alpha1[j]-alpha1[j-1])*180/3.14, (alpha2[j]-alpha2[j-1])*180/3.14, (alpha3[j]-alpha3[j-1])*180/3.14)
    j=j+1

j=0
while j<traj.n_frames:
    k=0
    while k<lag:
        if j+k<traj.n_frames:
            #     length1=[numpy.sum(axis11, e_vectors2[0]),numpy.dot(e_vectors1[1], e_vectors2[1]),numpy.dot(e_vectors1[2], e_vectors2[2])]
            #     cdot=numpy.dot(e_vectors1, numpy.transpose(e_vectors2))
            #     d[j+k]=traj.xyz[j+k,i2]-traj.xyz[j+k,i]
            #cdot=d[j][0]*d[j+k][0]+d[j][1]*d[j+k][1]+d[j][2]*d[j+k][2]
            #cos20=cdot1[0]/math.sqrt(cdot1[0]**2+cdot1[1]**2)
            #cos21=cdot2[0]/math.sqrt(cdot2[0]**2+cdot2[1]**2)
            #cos22=cdot3[0]/math.sqrt(cdot3[0]**2+cdot3[1]**2)
            corr0[k]=corr0[k]+(alpha1[j+k]-alpha1[j])**2   #math.acos(math.sqrt(cos20**2))**2 #0.5*(3*cos20**2-1)
            corr1[k]=corr1[k]+(alpha2[j+k]-alpha2[j])**2   #math.acos(math.sqrt(cos21**2))**2 #0.5*(3*cos21**2-1)
            corr2[k]=corr2[k]+(alpha3[j+k]-alpha3[j])**2   #math.acos(math.sqrt(cos22**2))**2 #0.5*(3*cos22**2-1)
            sum1[k]=sum1[k]+1
            #            print(k, alpha1[j]*180/3.14,(alpha1[j+k]-alpha1[j])*180/3.14)
        k=k+1
    j=j+1
RMASDfile=open(str(sys.argv[3]) + "RMASD.dat", 'w+')
#densityFILE1=open("../Data/correlation1opc.dat", 'w+')
#densityFILE2=open("../Data/correlation2opc.dat", 'w+')
k=0
while k<lag:
    print(k*timestep*0.001, corr0[k]/sum1[k], corr1[k]/sum1[k], corr2[k]/sum1[k],file=RMASDfile)
    k=k+1
    #    print(k, corr1[k]/sum1[k],file=densityFILE1)
    #    print(k, corr2[k]/sum1[k],file=densityFILE2)


#        i=i+1 
#        i2=i2+1
    
#hydrogens=traj.top.select("=~ 'H'")
#for i in hydrogens:
#    for  i2 in hydrogens:
        #i=29
#        atom_pairs=[[i , i2]]
#        r=md.compute_distances(traj,atom_pairs)
#        rAV=sum(r)/float(len(r))
#        if rAV<1.0 and i2>i:
#            j=0
#            while j<traj.n_frames:
#                d[j]=traj.xyz[j,i2]-traj.xyz[j,i]
#                k=0
                

