import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class image3d(object):
    '''
    Image3d
    Voxel scalar data
    '''
    pass

    def __init__(self,data,res=1):
        ''' 
		Build the image3d object
		:param data: array Nx3 dimention
		:type data: np.array
		:param res: resolution of the image (default 1)
		:type res: float
		'''
        self.im=data
        self.res=res
        return
		
    
    def xcorr3d(self,pad=1):
        '''
        3d autocorrelation
        :return: Autocor3d
        :rtype: image3d
        :return: Cinf
        :rtype: float
        '''
        
        import image3d.xcorr3d as xcorr3d
        mean_data=np.nanmean(self.im)
        
        if pad==1:
            An=np.fft.ifftn(np.abs(np.fft.fftn(self.im))**2)
            Autocor=np.abs(np.fft.fftshift(An/np.nanmax(An)))
            Cinf=mean_data**2/np.mean(self.im**2)
        elif pad==2:
            ss=np.shape(self.im)
            sn=np.array([ss[0]*pad,ss[1]*pad,ss[2]*pad])
            mpad=np.ones(sn)*mean_data
            mpad[0:ss[0],0:ss[1],0:ss[2]]=self.im
            An=np.fft.ifftn(np.abs(np.fft.fftn(mpad))**2)
            Autocor=np.abs(np.fft.fftshift(An/np.nanmax(An)))
            Cinf=mean_data**2/np.mean(mpad**2)
        else:
            return 'you should not do pad higher than 2, but if you think your machine can handle it do it yourself '
            
        return xcorr3d.xcorr3d(Autocor,self.res,Cinf)

		
	
    def extract_profil(self,init,end):
        '''
        Extract data along a line starting at init [xi,yi,zi] and ending at end [xe,ye,ze]
        '''

        radius=((end[0]-init[0])**2+(end[1]-init[1])**2+(end[2]-init[2])**2)**0.5
        vector=np.array([end[0]-init[0],end[1]-init[1],end[2]-init[2]])
        vector=vector/np.linalg.norm(vector)
		
        res=np.zeros(int(radius))
        xl=np.zeros(int(radius))
        for i in list(range(int(radius))):
            xi=init[0]+i*vector[0]
            yi=init[1]+i*vector[1]
            zi=init[2]+i*vector[2]
            xl[i]=((xi-init[0])**2+(yi-init[1])**2+(zi-init[2])**2)**0.5*self.res
            res[i]=self.im[int(xi),int(yi),int(zi)]

        return res,xl

    def plot(self,axis,pc,colorbar=cm.viridis):
	        '''
	        plot a slice of the 3d image
        
	        :param axis: slice perpendiculaire to the axis 'X','Y','Z'
	        :type axis: str
	        :param pc: ratio of the position between 0 to 1
	        :type pc: str
	        '''
        
	        ss=np.shape(self.im)
        
	        if (axis=='Z'):
        		img=plt.imshow(self.im[int(pc*(ss[0]-1)),:,:],cmap=colorbar,extent=(0,ss[0]*self.res,0,ss[1]*self.res))
        		plt.xlabel('+X')
        		plt.ylabel('-Y')
        	elif (axis=='Y'):
            		img=plt.imshow(self.im[:,int(pc*(ss[1]-1)),:],cmap=colorbar,extent=(0,ss[0]*self.res,0,ss[2]*self.res))
            		plt.xlabel('+X')
            		plt.ylabel('-Z')
        	elif (axis=='X'):
            		img=plt.imshow(self.im[:,:,int(pc*(ss[2]-1))],cmap=colorbar,extent=(0,ss[1]*self.res,0,ss[2]*self.res))
            		plt.xlabel('+Y')
            		plt.ylabel('-Z')
        
        	plt.axis('equal')
        	plt.colorbar(img,orientation='vertical',aspect=4)
        	return img