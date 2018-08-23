


# # CTMaging for Gate Simulation 

# ## Reading Dicom through VTK




import dicom
import os, sys
import numpy
import vtk
from vtk.util import numpy_support
from matplotlib import pyplot, cm
#get_ipython().magic(u'pylab inline')
#from natsort import natsorted
import plotly.offline as py
import plotly.graph_objs as go
py.offline.init_notebook_mode(connected=True)
from vtk.util import vtkImageImportFromArray as viifa 
import array


# ### Set the directory path of Dicom CT scann

# In[22]:

CTscan=sys.argv[1]

print
print "                 *****  CT scan Processing ****"
print 
print 
print "     * Dicom is being read ..."

PathDicom = CTscan
reader = vtk.vtkDICOMImageReader()
reader.SetDirectoryName(PathDicom)
reader.Update()


# ### Transform VTK object (CT scan into numpy 3d array)




#Get the 'vtkImageData' object from the reader
imageData = reader.GetOutput()

#Get the 'vtkPointData' forom the 'vtkImagedata' object
pointData = imageData.GetPointData()

#Ensure that only one array exist within the 'vtkPoint Data object
assert (pointData.GetNumberOfArrays()==1)

#get the 'vtkArray' (or whatever drived type) which is nedded for  the numpy_support.vtk_to_numpy` function
arrayData=pointData.GetArray(0)

#get the 'vtkArray' to Numpy array
ArrayDicom = numpy_support.vtk_to_numpy(arrayData)

#Get image dimension
DicomSize=reader.GetOutput().GetDimensions()

#Get pixel dimension
DicompixelSize=reader.GetDataSpacing()

#reshape the Numpy array 3D usin 'DicomSize' as a 'shape'
ArrayDicom = ArrayDicom.reshape(DicomSize, order='F')


print
print
print "     * Dicom Information ..."
print
print '            ** Dicom CT scan Dimensions : ' +str(DicomSize)
print '            ** Pixels Dimensions : ' + str(DicompixelSize)

# FLIP ARRAY TO HAVE A DORSAL RECUMBENCY 
arrayDicom=numpy.fliplr(ArrayDicom)



# ### Write a Raw File of CT scan 

# #### Image can be opened with image : 16 bit signed - Little endian byte order

# In[24]:


#data = np.zeros((DicomSize), dtype=np.int)
data = arrayDicom.flatten('F')
#print type(data)

print
print "     * Writing a Raw file of  Dicom CT scan..."
print


#Write Raw File 
filename = CTscan
binfile = filename + '.raw'
bitas = array.array('h', data)
#print binfile 
with open(binfile, 'wb') as file:
    file.write(bitas)

#Write mhd File
mhdFile = filename + '.mhd'
#print mhdFile

with open(mhdFile, 'wt') as file:
    file.write('ObjectType = Image\n')
    file.write('NDims = ' + str(len(DicomSize)) + '\n')
    file.write('BinaryData = True\n')
    file.write('BinaryDataByteOrderMSB = False\n')
    file.write('CompressedData = False\n')
    file.write('TransformMatrix = 1 0 0 0 1 0 0 0 1\n')
    file.write('Offset = 0 0 0\n')
    file.write('CenterOfRotation = 0 0 0\n')
    file.write('ElementSpacing = ' + str(DicompixelSize[0]) +' '+ str(DicompixelSize[1]) +' '+ str(DicompixelSize[2]) +'\n')
    file.write('DimSize = ' + str(DicomSize[0]) +' '+ str(DicomSize[1]) +' '+ str(DicomSize[2]) + '\n')
    file.write('ElementType = MET_SHORT\n')
    file.write('ElementDataFile = ' + binfile)
#print DicomSize

print '            ** Dicom .Raw file name : ' +str(binfile)
print '            ** Dicom .Mhd file name : ' +str(mhdFile)
print '            ** Raw file Dimensions : ' +str(DicomSize)


# ### CROP CT scan for GATE simulation 

# In[25]:

print
print "     * Cropping CT scan for GATE simulation..."
print
print "     * Writing a Raw file of  shrank CT scan..."
print

# Get index for cropping based on CT scan

#Create a sum of each slice of CT scan
SuM_CT=numpy.zeros((DicomSize[0], DicomSize[1]))
for k in range(0, DicomSize[2]):
    SuM_CT = SuM_CT + arrayDicom[: ,: , k]

#create meanvalue for Row (512, meanA) and for Column (512, meanB) 
meanA=numpy.mean(SuM_CT, axis=0)
meanB=numpy.mean(SuM_CT, axis=1)


#Define a threshlod 35% of the max HU for y index 56% of the HU for x index
thresholdA = max(meanA)/0.35
thresholdB = max(meanB)/0.56

#See curves bellow to set an addequat HU threshold for croping ( approx. 20% of max)
threshold = -125000

for i in range(0, len(meanA)):
    if meanA[i]>thresholdA:
        y1_crop=i-20  #ventre
        break

for i in range(len(meanA)-1, -1, -1):
    if meanA[i]>thresholdA:
        y2_croptemp=i-60  # suppression d'une partie de la table 
        break

for i in range(0, len(meanB)):
    if meanB[i]>thresholdB:
        x1_crop=i-20  #gauche de l'image
        break
        
for i in range(len(meanB)-1, -1, -1):
    if meanB[i]>thresholdB:
        x2_croptemp=i+20  #droite de l'image
        break
        

x1 = numpy.arange(0, len(meanA))
x2 =  numpy.arange(0, len(meanB))

#plot with plotly

#trace1 = go.Scatter(
#    x = x1,
#    y = meanA,
#    mode = 'lines', 
#    name = 'meanA')

#trace2 = go.Scatter(
#    x = x2,
#    y = meanB,
#    mode = 'lines', 
#    name = 'meanB')

#data = [trace1, trace2]
#py.iplot(data, filename='scatter-mode')

x2_crop=x2_croptemp-DicomSize[0]
y2_crop=y2_croptemp-DicomSize[1]


# crope numpy array of CT scann with array slicing [du iéme pixelà partir de la gauche:au iéme pixel à partir de ladroit, du iéme pixel à partir du haut:au iéme pixel à partir du bas]
crop_img=arrayDicom[x1_crop:x2_crop, y1_crop:y2_crop]
#print crop_img.shape

#Write raw file of copped CT scann 

datacrop = crop_img.flatten('F')

filename = CTscan+'_crop'
binfile = filename + '.raw'
bitas = array.array('h', datacrop)
#print binfile 
with open(binfile, 'wb') as file:
    file.write(bitas)

mhdFile = filename + '.mhd'
#print mhdFile

with open(mhdFile, 'wt') as file:
    file.write('ObjectType = Image\n')
    file.write('NDims = ' + str(len(crop_img.shape)) + '\n')
    file.write('BinaryData = True\n')
    file.write('BinaryDataByteOrderMSB = False\n')
    file.write('CompressedData = False\n')
    file.write('TransformMatrix = 1 0 0 0 1 0 0 0 1\n')
    file.write('Offset = 0 0 0\n')
    file.write('CenterOfRotation = 0 0 0\n')
    file.write('ElementSpacing = ' + str(DicompixelSize[0]) +' '+ str(DicompixelSize[1]) +' '+ str(DicompixelSize[2]) +'\n')
    file.write('DimSize = ' + str(crop_img.shape[0]) +' '+ str(crop_img.shape[1]) +' '+ str(crop_img.shape[2]) + '\n')
    file.write('ElementType = MET_SHORT\n')
    file.write('ElementDataFile = ' + binfile)
#print crop_img.shape

print '            ** Dicom .Raw file name : ' +str(binfile)
print '            ** Dicom .Mhd file name : ' +str(mhdFile)
print '            ** Raw file Dimensions : ' +str(crop_img.shape)
print


print "Raw images can be opened with ImageJ : 16 bits signed - Litlle endian byte order"
print "End !"
print
print "         * Jeremy Leste, Luc Simon, CRCT, janvier 2018"
print
print
