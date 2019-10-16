import struct
import numpy as np
import os
import math


#####################################################################################
#################################### INPUT:BEGIN ####################################

# path to the raw file
filename="somefile.raw"

# control whether to write the .perseus file
write_perseus_file=True

# control whether to write the raw file (necessary if the volume is downsampled)
write_vol_file=False

# the resolution of the volume
x_res=128
y_res=128
z_res=128

# control wether to downsample the volume
x_step=1
y_step=1
z_step=1

# endian of the data
endian="little"
# endian="big"

# data type of the voxels
datatype="uchar"
# datatype="float"
# datatype="double"

# number of bytes of each voxel
datasize=1

# whether to reverse the value of the voxels
reverse_val=False

# the following inputs should be ignored

case_val=None

storeorder="x_first"

#################################### INPUT:END ####################################
###################################################################################


x_cnt=int(x_res/x_step)
y_cnt=int(y_res/y_step)
z_cnt=int(z_res/z_step)

if case_val is not None:
    x_cnt+=2
    y_cnt+=2
    z_cnt+=2

if endian=="little":
    format="<"
elif endian=="big":
    format=">"

if datatype=="float":
    format+="f"
elif datatype=="double":
    format+="d"
elif datatype=="uchar":
    format+="B"

print("format: "+format)

# get offset of the point in the data array
def getOffset(x, y, z):
    if storeorder == "x_first":
        return z*z_step * x_res*y_res + y*y_step * x_res + x*x_step;
    elif storeorder == "z_first":
        return x*x_step * y_res*z_res + y*y_step * z_res + z*z_step;

with open(filename, mode='rb') as fin:
    filecontent = fin.read()

# print(type(filecontent))
# print(filecontent)

basename=os.path.basename(filename)
pos=basename.rfind('.')
if pos != -1:
    basename=basename[0:pos]

basename = basename + "_"+str(x_cnt)+"_"+str(y_cnt)+"_"+str(z_cnt)
if reverse_val:
    basename += "_rev"
if case_val is not None:
    basename += "_wcase"
outfname = basename + ".perseus"
outfname_vol = basename + ".raw"
# print(outfname)
# print(outfname_vol)

if write_perseus_file:
    fout=open(outfname, mode='w')
    fout.write("3\n");
    fout.write(str(x_cnt)+"\n")
    fout.write(str(y_cnt)+"\n")
    fout.write(str(z_cnt)+"\n")

if write_vol_file:
    fout_vol=open(outfname_vol, mode='wb')

for z in range(0,z_cnt):
    for y in range(0,y_cnt):
        for x in range(0,x_cnt):
            # print(str(x)+","+str(y)+","+str(z)+"\n")

            ind=None
            if case_val is not None:
                if x==0 or x==x_cnt-1 or y==0 or y==y_cnt-1 or z==0 or z==z_cnt-1:
                    scalar=case_val
                else:
                    ind=getOffset(x-1, y-1, z-1)
            else:
                ind=getOffset(x, y, z)

            if ind is not None:
                val=struct.unpack(format, filecontent[ind*datasize:(ind+1)*datasize])
                scalar=val[0]
                # scalar=math.log(val[0])
                if reverse_val:
                    scalar = -scalar

                if write_vol_file:
                    barr=struct.pack(format, scalar)
                    fout_vol.write(barr)

            if write_perseus_file:
                fout.write(str(scalar)+"\n")










