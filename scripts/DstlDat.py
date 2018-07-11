"""
    DstlDat.py

    Class to read in DSTL dat format files.

    Crown Copyright 2010

    For government use only.
    
"""
import numpy
import struct
import os

class DstlDat():
    def __init__(self, fname, mode = 'r'):
        self.fname = fname
        self.MAGIC = 0x87654321
        self.end_type = "<"
        
        self.int64_str = "Q"
        if struct.calcsize(self.int64_str) != 8:
            self.int64_str = "L"
            if struct.calcsize(self.int64_str) != 8:
                print "Can't figure out what an uint64_t should be read as!"
        if mode == 'r':
            fp = open(fname, "rb")
            magic_text =fp.read(4)

            magic = struct.unpack(self.end_type+"I", magic_text)[0]
            print "Magic = {0:x}".format(magic)
            if magic != self.MAGIC:
                self.end_type = ">"
                magic = struct.unpack(self.end_type+"I", magic_text)[0]
                print "Magic = {0:x}".format(magic)
                if magic != self.MAGIC:
                    print 
                    print "Failed to read in file!!!!"
                    return

            self.nx = struct.unpack(self.end_type+self.int64_str,fp.read(8))[0]
            self.ny = struct.unpack(self.end_type+self.int64_str,fp.read(8))[0]

            self.sx = struct.unpack(self.end_type+"d",fp.read(8))[0]
            self.sy = struct.unpack(self.end_type+"d",fp.read(8))[0]

            self.itype = struct.unpack(self.end_type+"i",fp.read(4))[0]
            print str(self.nx) + " X " + str(self.ny)

            print "data type = " + str(self.itype)
            if self.itype == 2:
                self.dtype = numpy.complex64
                self.dtype_a = "f"
                self.dsize = 8

            if self.itype == 3:
                self.dtype = numpy.float
                self.dtype_a = "f"
                self.dsize = 4

            if self.itype == 4:
                self.dtype = numpy.double
                self.dtype_a = "d"
                self.dsize = 8

            if self.itype == 6:
                self.dtype = numpy.int32
                self.dtype_a = "i"
                self.dsize = 4
            
            if self.itype == 13:
                self.dtype = numpy.double
                self.dtype_a = 'd'
                self.dsize = 3 * 8

            if self.itype == 11:
                self.dtype = numpy.int16
                self.dtype_a = 'h'
                self.dsize = 2
                
        else:
            fp = open(fname, "wb")
            
        fp.close()
            
    def load(self, ibnds=(-1,-1,-1,-1)):
        """ Load up data from the DstlDat file. Pass in the bounds as a tuple:
            (start_row, end_row, start_col, end_col)
            """
        
        fp = open(self.fname, "rb")

        if ibnds[0] == -1 and ibnds[1] == -1 and ibnds[2] == -1 and ibnds[3] == -1:            
            bnds = (0, self.ny, 0, self.nx)
        else:
            bnds = ibnds

        
        nr = bnds[1]-bnds[0]
        nc = bnds[3]-bnds[2]

        if nr > self.ny or bnds[1] > self.ny:
            raise IndexError('Number of rows requested bigger than that stored in file!')

        if nc > self.nx or bnds[3] > self.nx:
            raise IndexError('Number of columns requested bigger than that stored in file!')
        
        if self.itype == 13:
            if nr == 1 or nc == 1:
                data = numpy.zeros((nr*nc,3), dtype=self.dtype)
            else:
                data = numpy.zeros((nr,nc,3), dtype=self.dtype)
        else:
            if nr == 1 or nc == 1:
                data = numpy.zeros((nr*nc,), dtype=self.dtype)
            else:
                data = numpy.zeros((nr,nc), dtype=self.dtype)
        
        z = 1j
        fp.seek(40 + ((bnds[0] * self.nx + bnds[2]) * self.dsize ))
        
        if self.dtype == numpy.complex64:
            for rows in range(nr):
                ldat = struct.unpack(self.end_type+str(nc*2)+self.dtype_a, fp.read(self.dsize * nc))
                if nr == 1:
                    data[:] = numpy.array(ldat[::2]) + z * numpy.array(ldat[1::2])
                elif nc == 1:
                    data[rows] = ldat[0] + z * ldat[1]
                else:
                    data[rows,:] = numpy.array(ldat[::2]) + z * numpy.array(ldat[1::2])
                fp.seek(((self.nx - bnds[3]) + bnds[2]) * self.dsize, os.SEEK_CUR)
        elif self.itype == 13:
            for rows in range(nr):
                ldat = numpy.array(struct.unpack(self.end_type+str(3*nc)+self.dtype_a, fp.read(self.dsize * nc)))
                if nr == 1:
                    data =  numpy.reshape(ldat, data.shape)
                elif nc == 1:
                    data[rows] = numpy.reshape(ldat, (nc, 3))
                else:
                    data[rows,:, :] =numpy.reshape(ldat, (nc, 3))
                fp.seek(((self.nx - bnds[3]) + bnds[2]) * self.dsize, os.SEEK_CUR)
        else:
            for rows in range(nr):
                ldat = fp.read(self.dsize * nc)
                if nr == 1 or nc == 1:
                    data[:] = struct.unpack(self.end_type+str(nc)+self.dtype_a, ldat)
                elif nc == 1:
                    data[rows] = struct.unpack(self.end_type+str(nc)+self.dtype_a, ldat)[0]
                else:
                    data[rows,:] = struct.unpack(self.end_type+str(nc)+self.dtype_a, ldat)
                fp.seek(((self.nx - bnds[3]) + bnds[2]) * self.dsize, os.SEEK_CUR)
        fp.close()

        return (data)
        
    def save(self, data, parameter): # parameter is dictionary containing fields sx and sy
        fp = open(self.fname, "w+b")
        magic = struct.pack(self.end_type+"I", self.MAGIC)
        parameter["nx"] = data.shape[1]
        parameter["ny"] = data.shape[0]
        fp.write(magic) # Write 4 bytes

        fp.write(struct.pack(self.end_type+self.int64_str, parameter["nx"])) # Write 8 bytes
        fp.write(struct.pack(self.end_type+self.int64_str, parameter["ny"])) # Write 8 bytes
        
        fp.write(struct.pack(self.end_type+"d", parameter["sx"])) # Write 8 bytes
        fp.write(struct.pack(self.end_type+"d", parameter["sy"])) # Write 8 bytes
        
        dataType = data.dtype
        if dataType == "complex64":
            self.dtype = 2
            self.dtype_a = "f"
            self.dsize = 8
        elif data.shape[-1] == 3 and dataType == "float64":
            self.dtype = 13
            self.dtype_a = 'd'
            self.dsize = 3 * 8
        elif dataType == "float32":
            self.dtype = 3
            self.dtype_a = "f"
            self.dsize = 4
        elif dataType == "float64":
            self.dtype = 4
            self.dtype_a = "d"
            self.dsize = 8
        elif dataType == "int32":
            self.dtype = 6
            self.dtype_a = "i"
            self.dsize = 4
        else:
            print "data type " + dataType + " not supported."
            return
        fp.write(struct.pack(self.end_type+"i", self.dtype)) # Write 4 bytes
        
        # Now 40 bytes into file
        dataRow = numpy.array(parameter["nx"]*2*[1], dtype = dataType)
        if dataType == "complex64":
            for rows in range(parameter["ny"]):
                dataRow[::2] = data[rows, :].real
                dataRow[1::2] = data[rows, :].imag
                fp.write(dataRow.tostring())
        elif self.dtype == 13:
            for rows in range(parameter["ny"]):
                dataRow = data[rows, :, :].reshape(parameter["nx"]*3)
                fp.write(dataRow.tostring())
        else:
            for rows in range(parameter["ny"]):
                fp.write(data[rows, :].tostring())
        fp.close()
        
if __name__ == "__main__":    
    data = numpy.array([[1.3, 2.4, 5.6] * 4, [1.2, 3.5, 0.7, 0.4]*3], dtype = "float64")
    data = data.reshape([2, 4, 3])
    param = {"sx" : 0.4, "sy" : 0.6}
    filePoint = DstlDat("test.dat", "w")
    filePoint.save(data, param)
    
    filePoint2 = DstlDat("test.dat")
    data2 = filePoint2.load()
    if (data-data2).any() != 1:
        print "Testing Completed Successfully"
    os.remove("test.dat")


        
