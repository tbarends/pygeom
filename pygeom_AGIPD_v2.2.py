
################################################ USAGE


# Use these statements on psexport
#infilename=sys.argv[1]
#geomfilename=sys.argv[2]
#maxpredictions=int(sys.argv[3])

# Use these statements on Thomas' laptop for testing
#infilename='lyso-gd-4.stream'
#geomfilename='stephan-kn4.geom'
#workingdir='C:/cygwin64/home/tbarend'
#os.chdir(workingdir)
#maxpredictions=300000


################################################

#import modules
import numpy as np
import scipy.stats as st
import matplotlib
matplotlib.use("tkagg")
import matplotlib.pyplot as plt
from matplotlib.pylab import *
import string
import sys

########################READGEOMFILE#######################################
def readgeomfile(infilename):
    infile=open(infilename,'r')
    panels=[]
    clen=''
    offset=''
    adu=''
    photonenergy=''
    geometry=np.zeros((128,10))
    panelname='humbugtotalhumbugnobodywillcallapanellikethis'
    panel=-1
    for line in infile:
        truncated=line[0:len(line)-1]
        elements=truncated.split('=')
        if 'clen' in line:
            clen=line
            continue
        if 'offset' in line:
            offset=line
            continue
        if 'adu' in line:
            adu=line
            continue
	if 'photon_energy' in line:
            photonenergy=line
            continue
	    
        if not '/' in line or ";" in line:
            continue
       



        
        identifiers=elements[0].split('/')
        

        if identifiers[0]<>panelname:
            panel=panel+1
            #print 'found panel',identifiers[0]
            panelname=identifiers[0]
            panels.append(panelname)
        if 'min_fs' in identifiers[1]:
            geometry[panel,0]=int(float(elements[1]))
            continue
        if 'min_ss' in identifiers[1]:
            geometry[panel,1]=int(float(elements[1]))
            continue
        if 'max_fs' in identifiers[1]:
            geometry[panel,2]=int(float(elements[1]))
            continue
        if 'max_ss' in identifiers[1]:
            geometry[panel,3]=int(float(elements[1]))
            continue
        
        if 'corner_x' in identifiers[1]:
            geometry[panel,4]=(float(elements[1]))
            continue
        if 'corner_y' in identifiers[1]:
            geometry[panel,5]=(float(elements[1]))
            continue

        if 'fs' in identifiers[1]:
            functionleft=elements[1].split('x')
            geometry[panel,6]=float(functionleft[0])
            functionright=functionleft[1].split('y')
            geometry[panel,7]=float(functionright[0])

        if 'ss' in identifiers[1]:
            functionleft=elements[1].split('x')
            geometry[panel,8]=float(functionleft[0])
            functionright=functionleft[1].split('y')
            geometry[panel,9]=float(functionright[0])
    

    infile.close()
    return geometry,panels,adu,clen,offset,photonenergy

##########################WRITEGEOMFILE####################################
def writegeomfile(outfilename,geom,panels,adu,clen,offset,photonenergy):
    outfile=open(outfilename,'w')
    outfile.write(adu+'\n')
    outfile.write(clen+'\n')
    outfile.write(offset+'\n')
    outfile.write(photonenergy+'\n')
    outfile.write(' \n')
                  

        
    s=np.shape(geom)
    numpanels=s[0]
    for n in range(numpanels):
        identifier=panels[n]+'/'
        outfile.write(identifier+'min_fs = '+str(int(geom[n,0]))+'\n')
        outfile.write(identifier+'min_ss = '+str(int(geom[n,1]))+'\n')
        outfile.write(identifier+'max_fs = '+str(int(geom[n,2]))+'\n')
        outfile.write(identifier+'max_ss = '+str(int(geom[n,3]))+'\n')
        outfile.write(identifier+'badrow_direction = - \n')
        

        #fsxfactor=str(geom[n,6]).format('d')
        #fsyfactor=str(geom[n,7]).format('d')
        fsxfactor="%1.7f" % geom[n,6]
        fsyfactor="%1.7f" % geom[n,7]   
        fsformula=fsxfactor+'x'+' +'+fsyfactor+'y'
        if (geom[n,7])<0:
            fsformula=fsxfactor+'x'+' '+fsyfactor+'y'

        #ssxfactor=str(geom[n,8]).format('d')
        #ssyfactor=str(geom[n,9]).format('d')
        ssxfactor="%1.7f" % geom[n,8]                    
        ssyfactor="%1.7f" % geom[n,9]                                       
   
        ssformula=ssxfactor+'x'+' +'+ssyfactor+'y'
        if (geom[n,9])<0:
            ssformula=ssxfactor+'x'+' '+ssyfactor+'y'

        outfile.write(identifier+'fs = '+fsformula+'\n')
        outfile.write(identifier+'ss = '+ssformula+'\n')

        outfile.write(identifier+'corner_x = '+str(geom[n,4])+'\n')
        outfile.write(identifier+'corner_y = '+str(geom[n,5])+'\n')
       

        outfile.write(identifier+'no_index = 0 \n')
        outfile.write(' \n')
    outfile.close





########################MAPGEOMETRY########################################
def mapgeometry(geometry):
    numpanels=int(np.shape(geometry)[0]) # how many panels are there
    minfs_overall=int(np.min(geometry[:,0]))  # min/max indices for detector
    minss_overall=int(np.min(geometry[:,1]))
    maxfs_overall=int(np.max(geometry[:,2]))
    maxss_overall=int(np.max(geometry[:,3]))

    mapx=np.zeros((maxfs_overall+1,maxss_overall+1)) # initialize map
    mapy=np.zeros((maxfs_overall+1,maxss_overall+1)) # initialize map

    
    for panel in range(numpanels):       # loop over panels
        minfs=int(geometry[panel,0])
        minss=int(geometry[panel,1])
        maxfs=int(geometry[panel,2])
        maxss=int(geometry[panel,3])
        cornerx=float(geometry[panel,4])
        cornery=float(geometry[panel,5])
        fsxfactor=float(geometry[panel,6])
        fsyfactor=float(geometry[panel,7])
        ssxfactor=float(geometry[panel,8])
        ssyfactor=float(geometry[panel,9])
        
        fs=np.arange(0,maxfs-minfs+1)
        ss=np.arange(0,maxss-minss+1)
        
        panelpixels=np.meshgrid(fs,ss)
        #panelpixels contains 2 matrices! 0 has the fs mesh, 1 has the ss mesh

        x=cornerx+(fsxfactor*panelpixels[0])+(ssxfactor*panelpixels[1])
        y=cornery+(fsyfactor*panelpixels[0])+(ssyfactor*panelpixels[1])
        mapx[minfs:maxfs+1,minss:maxss+1]=np.transpose(x)
        mapy[minfs:maxfs+1,minss:maxss+1]=np.transpose(y)
        
    #print x
    #print y
    print 'Prepared pixelmaps of shapes',np.shape(mapx),np.shape(mapy)
    return mapx,mapy
    
########################PIXEL2LAB##############################################

def pixel2lab(fs,ss,geometry):
    
    numpanels=int(np.shape(geometry)[0]) # how many panels are there
    
    for panel in range(numpanels):       # loop over panels
        minfs=int(geometry[panel,0])
        minss=int(geometry[panel,1])
        maxfs=int(geometry[panel,2])
        maxss=int(geometry[panel,3])
        if (fs<=maxfs) and (fs>=minfs) and (ss<=maxss) and (ss>=minss):
            cornerx=float(geometry[panel,4])
            cornery=float(geometry[panel,5])
            fsxfactor=float(geometry[panel,6])
            fsyfactor=float(geometry[panel,7])
            ssxfactor=float(geometry[panel,8])
            ssyfactor=float(geometry[panel,9])
            x=cornerx + (fs*fsxfactor) + (ss*ssxfactor)
            y=cornery + (fs*fsyfactor) + (ss*ssyfactor)
            break
    return x,y


###############################MAIN############################################

# Use these statements on psexport
infilename=sys.argv[1]
geomfilename=sys.argv[2]
maxpredictions=int(sys.argv[3])

# Use these statements on Thomas' laptop for testing
#infilename='lyso-gd-4.stream'
#geomfilename='stephan-kn4.geom'
#workingdir='C:/cygwin64/home/tbarend'
#os.chdir(workingdir)
#maxpredictions=300000

matchcounter=0
numindexed=0
maxpeaks = 4000
cutoff = 10
#declare numpy arrays
peakx=np.zeros((maxpeaks))
peaky=np.zeros((maxpeaks))
predx=np.zeros((maxpredictions))
predy=np.zeros((maxpredictions))
pixfs=np.zeros((maxpredictions))
pixss=np.zeros((maxpredictions))
distx=np.zeros((maxpredictions))
disty=np.zeros((maxpredictions))
dist=np.zeros((maxpredictions))

# Read the geometry file
geom,panels,adu,clen,offset,photonenergy=readgeomfile(geomfilename)

mapx,mapy=mapgeometry(geom)

# Open the stream file
infile=open(infilename,'r')
print 'Opened',infilename,'for reading...'

# Flags for use inside loop 
indexed=0
inpeaklist=0
inpredlist=0

# Counters for numbers of peaks and predictions encountered
numpeak=0
numpred=0

#iterate over lines in the file
for line in infile:
    #print indexed
    if ('CrystFEL' in line):
        continue 
 
    if ('indexed_by = mosflm' in line) or ('indexed_by = dirax'):
        #print line
        indexed=1
        numindexed=numindexed+1
        #continue

    if ('End chunk' in line):
        #print "End chunk"
        indexed=0
        numpeak=0
        continue
     
    if ('fs/px   ss/px (1/d)/nm^-1   Intensity  Panel' in line) and (indexed==1):
        inpeaklist=1
        #print "Found peaks for indexed image"
        continue

    if ('End of peak list' in line):
        #print "End of peak list encountered"
        inpeaklist=0
        continue

    if (inpeaklist==1):
        elements=line.split()
        #peakx[numpeak]=float(elements[0])
        #peaky[numpeak]=float(elements[1])
        #peakx[numpeak],peaky[numpeak]=pixel2lab(float(elements[0]),float(elements[1]),geom)
        peakx[numpeak]=mapx[int(float(elements[0])),int(float(elements[1]))]
        peaky[numpeak]=mapy[int(float(elements[0])),int(float(elements[1]))]


        
        numpeak=numpeak+1
        continue

    if ('h    k    l' in line) or (' h   k   l' in line):
        inpredlist=1
        #print "Found a prediction list!"
        continue
    
    if ('End of reflections' in line):
        inpredlist=0
        indexed=0
        continue

    if (inpredlist==1):
        elements=line.split()
        #x=float(elements[7])
        #y=float(elements[8])
        #x,y=pixel2lab(float(elements[7]),float(elements[8]),geom)

        x=mapx[int(float(elements[7])),int(float(elements[8]))]
        y=mapy[int(float(elements[7])),int(float(elements[8]))]
               
        
        #predx[numpred]=x
        #predy[numpred]=y
        #pixx[numpred]=float(elements[7])
        #pixy[numpred]=float(elements[8])
        #print x,y
        bestdist=1000000
        bestpeakx=0
        bestpeaky=0
        for n in range(numpeak):        
            distance=np.sqrt(np.square(x-peakx[n])+np.square(y-peaky[n]))
            #print distance
            if (distance<bestdist):
                bestdist=distance
                bestpeakx=peakx[n]
                bestpeaky=peaky[n]
        if (bestdist<cutoff):
            #print "Found a prediction at",x,y,"close to a peak at",bestpeakx,bestpeaky,'at a distance of',bestdist
            dist[numpred]=bestdist
            distx[numpred]=x-bestpeakx
            disty[numpred]=y-bestpeaky

            predx[numpred]=x
            predy[numpred]=y
            pixfs[numpred]=float(elements[7])
            pixss[numpred]=float(elements[8])
            matchcounter=matchcounter+1
            if matchcounter==100:
                print "Found",numpred+1,"matches of a target of",maxpredictions,"(",float(100.0*(numpred+1.0)/maxpredictions),"% done)"
                matchcounter=0
            numpred=numpred+1
            if (numpred>maxpredictions-1):
                break

        



print "Found",numpred,"predictions in file"
print "Found",numpeak,"peaks in last indexed image"
print "The distances are between",np.min(dist),"and",np.max(dist)
print "The average discrepancy in x over the whole detector is",np.mean(distx)
print "The average discrepancy in y over the whole detector is",np.mean(disty)
print "Discrepancies are calculated as X_prediction-X_peak or Y_prediction-Y_peak."
print "So the discrepancies should be added to the corner positions"
print "The predictions are between",np.min(predx),np.max(predx),"and ",np.min(predy),np.max(predy)
print " "
print "Now looking at the geometry file:"
s= np.shape(geom)
newgeom=geom
print "There appear to be",s[0],"ASICs"
for asic in range(s[0]):
    minfs=geom[asic,0]
    minss=geom[asic,1]
    maxfs=geom[asic,2]
    maxss=geom[asic,3]
    cornerx=geom[asic,4]
    cornery=geom[asic,5]
    print "Now working on ASIC #",asic
    #print "fs between",minfs,"and",maxfs
    #print "ss between",minss,"and",maxss
    #print "Corner x,y",cornerx,",",cornery
    #print "Going through the list of matched predictions:"
    xdiscr=[]
    ydiscr=[]
    x=[]
    y=[]
    fs=[]
    ss=[]
    
    for n in range(numpred):
        if (pixfs[n]>minfs and pixfs[n]<maxfs and pixss[n]>minss and pixss[n]<maxss):
            #print "Found a prediction for this asic at fs,ss", pixfs[n],pixss[n],"with x,y discrepancy",distx[n],disty[n]
            xdiscr.append(distx[n])
            ydiscr.append(disty[n])
            x.append(predx[n])
            y.append(predy[n])
            fs.append(pixfs[n])
            ss.append(pixss[n])

            
    if len(xdiscr)==0:
        print "No matches found on this ASIC - skipping"
        print "   "
        continue
    print "Found", len(xdiscr), 'matches'   
    print "Average discrepancies in x,y for this asic",np.average(np.asarray(xdiscr)),np.average(np.asarray(ydiscr))
    print "Median discrepcancies in x,y",np.median(np.asarray(xdiscr)),np.median(np.asarray(ydiscr))
    print "Standard deviations of discrepancies      ",np.std(np.asarray(xdiscr)),np.std(np.asarray(ydiscr))
    newgeom[asic,4]=cornerx+np.average(np.asarray(xdiscr))
    newgeom[asic,5]=cornery+np.average(np.asarray(ydiscr))

    if len(xdiscr)<50:
        print "Too few matches on this ASIC - skipping determination of rotations"
        print "   "
        continue

    fsxc=np.polyfit(fs,xdiscr,1)
    fsxcorr=st.pearsonr(fs,xdiscr)
    print "The dependence of the x discrepancy on fs is dx=",fsxc[0],'.fs+',fsxc[1],', correlation=',fsxcorr[0]
       
    fsyc=np.polyfit(fs,ydiscr,1)
    fsycorr=st.pearsonr(fs,ydiscr)
    print "The dependence of the y discrepancy on fs is dy=",fsyc[0],'.fs+',fsyc[1],', correlation=',fsycorr[0]

    ssxc=np.polyfit(ss,xdiscr,1)
    ssxcorr=st.pearsonr(ss,xdiscr)
    print "The dependence of the x discrepancy on ss is dx=",ssxc[0],'.ss+',ssxc[1],', correlation=',ssxcorr[0]

    ssyc=np.polyfit(ss,ydiscr,1)
    ssycorr=st.pearsonr(ss,ydiscr)
    print "The dependence of the y discrepancy on ss is dy=",ssyc[0],'.s+',ssyc[1],', correlation=',ssycorr[0]

    #if abs(fsxcorr[0])>0.2:
    newgeom[asic,6]=newgeom[asic,6]+(0.3333333*fsxc[0])
        
    #if abs(fsycorr[0])>0.2:
    newgeom[asic,7]=newgeom[asic,7]+(0.3333333*fsyc[0])

    #if abs(ssxcorr[0])>0.2:
    newgeom[asic,8]=newgeom[asic,8]+(0.3333333*ssxc[0])

    #if abs(ssycorr[0])>0.2:
    newgeom[asic,9]=newgeom[asic,9]+(0.3333333*ssyc[0])
    

##ddd    xxc=np.polyfit(x,xdiscr,1)
##    xxcorr=st.pearsonr(x,xdiscr)
##    print "The dependence of the x discrepancy on x is dx=",xxc[0],'.x+',xxc[1],', correlation=',xxcorr[0]
##       
##    xyc=np.polyfit(x,ydiscr,1)
##    xycorr=st.pearsonr(x,ydiscr)
##    print "The dependence of the y discrepancy on x is dy=",xyc[0],'.x+',xyc[1],', correlation=',xycorr[0]
##
##    yxc=np.polyfit(y,xdiscr,1)
##    yxcorr=st.pearsonr(y,xdiscr)
##    print "The dependence of the x discrepancy on y is dx=",yxc[0],'.y+',yxc[1],', correlation=',yxcorr[0]
##
##    yyc=np.polyfit(y,ydiscr,1)
##    yycorr=st.pearsonr(y,ydiscr)
##    print "The dependence of the y discrepancy on y is dy=",yyc[0],'.y+',yyc[1],', correlation=',yycorr[0]

    print "   "
    
print 'Writing updated geometry file output.geom'
writegeomfile('output.geom',newgeom,panels,adu,clen,offset,photonenergy)
    



font = {'color'  : 'black',
        'weight' : 'normal',
        'size'   : 10,
        }

print "Writing plot files..."
fig,ax=plt.subplots()
cax=ax.scatter(predx[0:numpred],predy[0:numpred],c=dist[0:numpred])
cbar=fig.colorbar(cax)
plt.title('Distance between predictions and closest peak mapped on detector area',fontdict=font)
plt.savefig('totaldistance.png')

fig,ax=plt.subplots()
cax=ax.scatter(predx[0:numpred],predy[0:numpred],c=distx[0:numpred])
cbar=fig.colorbar(cax)
plt.title('Horizontal distance between predictions and closest peak mapped on detector area',fontdict=font)
plt.savefig('xdistance.png')

fig,ax=plt.subplots()
cax=ax.scatter(predx[0:numpred],predy[0:numpred],c=disty[0:numpred])
plt.title('Vertical distance between predictions and closest peak mapped on detector area',fontdict=font)
cbar=fig.colorbar(cax)
plt.savefig('ydistance.png')

print "Showing plot files..."

fig,ax=plt.subplots()
cax=ax.scatter(predx[0:numpred],predy[0:numpred],c=dist[0:numpred])
cbar=fig.colorbar(cax)
plt.title('Distance between predictions and closest peak mapped on detector area',fontdict=font)
plt.show()

fig,ax=plt.subplots()
cax=ax.scatter(predx[0:numpred],predy[0:numpred],c=distx[0:numpred])
cbar=fig.colorbar(cax)
plt.title('Horizontal distance between predictions and closest peak mapped on detector area',fontdict=font)
plt.show()

fig,ax=plt.subplots()
cax=ax.scatter(predx[0:numpred],predy[0:numpred],c=disty[0:numpred])
plt.title('Vertical distance between predictions and closest peak mapped on detector area',fontdict=font)
cbar=fig.colorbar(cax)
plt.show() 

##print "Now assembling three arrays"
##im=np.zeros((1480,1552))
##imx=np.zeros((1480,1552))
##imy=np.zeros((1480,1552))
##
##
##
##for n in range(numpred-1):
##    x=int(predx[n])
##    y=int(predy[n])
##    if (x<1553) and (y<1481):
##        im[y,x]=dist[n]
##        imx[y,x]=distx[n]
##        imy[y,x]=disty[n]
##
##fig2,ax2=plt.subplots()
##cax2=ax2.imshow(im)
##plt.show()
##
##fig2,ax2=plt.subplots()
##cax2=ax2.imshow(imx)
##plt.show()
##
##fig2,ax2=plt.subplots()
##cax2=ax2.imshow(imy)
##plt.show()
##
##
##
##print "Now writing file in hdf5 format"
##outfile=h5py.File('distances.h5')
##distdataset=outfile.create_dataset('distances',data=im[0:1480,0:1552])
##distxdataset=outfile.create_dataset('distancesx',data=imx[0:1480,0:1552])
##distydataset=outfile.create_dataset('distancesy',data=imy[0:1480,0:1552])
##
##print im[0:1480,0:1552] 

print "all done!"


         

        

