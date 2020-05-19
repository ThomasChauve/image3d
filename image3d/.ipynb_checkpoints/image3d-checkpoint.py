import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import skimage
import setvector3d.setvector3d as sv3d
import scipy
import pyfftw

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
    
    
    def xcorr3d(self,pad=1,rad_Tukey=0,gray_level=True,usePYFFTW=False,nb_core=2):
        '''
        3d autocorrelation
        
        :param pad: Pading value 1 or 2
        :type pad: int
        :param rad_Tukey: radius for the Tukey function (pixel) (attenuation function : https://fr.mathworks.com/help/signal/ref/tukeywin.html)
        :type rad_Tukey: int
        :return: Autocorrelation function
        :rtype: xcorr3d.xcorr3d
        '''
        pyfftw.config.NUM_THREADS=nb_core
        
        import image3d.xcorr3d as xcorr3d
        mean_data=np.nanmean(self.im)
        
        if rad_Tukey!=0:
            ss=np.shape(self.im)
            coeff=tukeywin3D(ss,rad_Tukey)
            map=self.im*coeff
        else:
            map=self.im
        
        if gray_level:
            if pad==1:
                mpad=map
            elif pad==2:
                ss=np.shape(map)
                sn=np.array([ss[0]*pad,ss[1]*pad,ss[2]*pad])
                mpad=np.ones(sn)*mean_data
                mpad[0:ss[0],0:ss[1],0:ss[2]]=self.im
            elif (pad>1)and(pad<2):
                ss=np.shape(map)
                sn=np.array([int(ss[0]*pad),int(ss[1]*pad),int(ss[2]*pad)])
                mpad=np.ones(sn)*mean_data
                mpad[0:ss[0],0:ss[1],0:ss[2]]=self.im
            else:
                return 'you should not do pad higher than 2, but if you think your machine can handle it do it yourself '
            
            if usePYFFTW:
                An=pyfftw.interfaces.numpy_fft.ifftn(np.abs(pyfftw.interfaces.numpy_fft.fftn(map))**2)
            else:
                An=np.fft.ifftn(np.abs(np.fft.fftn(map))**2)
                
            Autocor=np.abs(np.fft.fftshift(An/np.nanmax(An)))
            Cinf=mean_data**2/np.mean(map**2)
            
        return xcorr3d.xcorr3d(Autocor,self.res,Cinf)
    

		
	
    def extract_profil(self,init,end):
        '''
        Extract data along a line starting at init [xi,yi,zi] and ending at end [xe,ye,ze]
        
        :param init: array dimention 3
        :type init: np.array
        :param end: array dimention 3
        :type end: np.array
        :return: res value along the profil
        :rtype: np.array
        :return: xl position along the profil
        :rtype: np.array
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
        
	        if (axis=='X'):
        		img=plt.imshow(self.im[int(pc*(ss[0]-1)),:,:],cmap=colorbar,extent=(0,ss[2]*self.res,0,ss[1]*self.res),origin='lower')
        		plt.xlabel('Z')
        		plt.ylabel('Y')
        	elif (axis=='Y'):
            		img=plt.imshow(self.im[:,int(pc*(ss[1]-1)),:],cmap=colorbar,extent=(0,ss[2]*self.res,0,ss[0]*self.res),origin='lower')
            		plt.xlabel('Z')
            		plt.ylabel('X')
        	elif (axis=='Z'):
            		img=plt.imshow(self.im[:,:,int(pc*(ss[2]-1))],cmap=colorbar,extent=(0,ss[1]*self.res,0,ss[0]*self.res),origin='lower')
            		plt.xlabel('Y')
            		plt.ylabel('X')
        
        	plt.axis('equal')
        	plt.colorbar(img,orientation='vertical',aspect=4)
        	return img
        
    def inertia_tensor(self):
        '''
        Compute the inertia tensor of the image and retrun the eigenvalue and eigenvector.
        :return eigval: eigenvalue
        :return eigvec: eigenvector
        '''
        inertia=skimage.measure.inertia_tensor(self.im)
        eigval,eigvec=np.linalg.eig(inertia)
            
        return eigval/self.res,eigvec
    
    def split_img(self,size_box,option):
        '''
        Divide the image in sub-image of box size size_box³
        :param size_box: size of the cubic box
        :type size_box: np.float
        :param option: (1) paved surface without the excess surface and without overlap (2) paved and with overlap without excess surface
        :type option: int
        '''
        center_box=[]
        sub_img=[]
        ss=np.shape(self.im)
        pix_sb=int(size_box/self.res)
        nbbox=np.array(ss)/pix_sb
        if option==1:
            xl=np.int64(np.linspace(0,pix_sb*int(nbbox[0]),int(nbbox[0]+1)))
            yl=np.int64(np.linspace(0,pix_sb*int(nbbox[1]),int(nbbox[1]+1)))
            zl=np.int64(np.linspace(0,pix_sb*int(nbbox[2]),int(nbbox[2]+1)))       
        elif option==2:
            return print('not implemented yet'),print('TO DO')
            
        for i in list(range(len(xl[0:-1]))):
            for j in list(range(len(yl[0:-1]))):
                for k in list(range(len(zl[0:-1]))):
                    center_box.append(np.array([xl[i]+pix_sb/2,yl[j]+pix_sb/2,zl[k]+pix_sb/2]))
                    sub_img.append(image3d(self.im[xl[i]:xl[i+1]-1,yl[j]:yl[j+1]-1,zl[k]:zl[k+1]-1],self.res))

        return center_box,sub_img
    
    def texture_anisotropy(self,size_box):
        '''
        Compute the inertia tensor on all the image for sub image of size size_box
        '''
        list_dict=[]
        for sbb in size_box:
            center_box,sub_img=self.split_img(sbb,1)
            eigval=[]
            eigvec=[]
            mainvec=[]
            cbox=[]
            center_mass=[]
            mass=[]
            
            # Anisotropy image
            tmp=np.zeros([len(center_box),3])
            for i in list(range(len(center_box))):
                tmp[i,:]=center_box[i]

            idx=len(np.unique(tmp[:,0]))
            xx=np.unique(tmp[:,0])
            idy=len(np.unique(tmp[:,1]))
            yy=np.unique(tmp[:,0])
            idz=len(np.unique(tmp[:,2]))
            zz=xx=np.unique(tmp[:,0])
            
            Ani_im_RA=np.zeros([idx,idy,idz])
            Ani_im_FA=np.zeros([idx,idy,idz])
            Ani_im_VA=np.zeros([idx,idy,idz])
            Ani_im_FlA=np.zeros([idx,idy,idz])
            
            
            for i in list(range(len(center_box))):
                iix=np.where(xx==center_box[i][0])[0]
                iiy=np.where(yy==center_box[i][1])[0]
                iiz=np.where(zz==center_box[i][2])[0]
                if np.sum(sub_img[i].im[:])!=0:
                    tmp_eigval,tmp_eigvec=sub_img[i].inertia_tensor()
                    cbox.append(center_box[i]*self.res)
                    center_mass.append((center_box[i]/self.res-sbb/2+scipy.ndimage.measurements.center_of_mass(sub_img[i].im)*self.res))
                    mass.append(np.sum(sub_img[i].im[:]))
                    eigval.append(tmp_eigval)
                    eigvec.append(tmp_eigvec)
                    id=np.where(tmp_eigval==np.min(tmp_eigval))[0][0]
                    mainvec.append(tmp_eigvec[:,id])
                    Ani_im_RA[iix,iiy,iiz]=np.std(tmp_eigval)/np.mean(tmp_eigval)
                    Ani_im_FA[iix,iiy,iiz]=np.std(tmp_eigval)/np.mean(np.multiply(tmp_eigval,tmp_eigval))**0.5
                    Ani_im_VA[iix,iiy,iiz]=1.-tmp_eigval[0]*tmp_eigval[1]*tmp_eigval[2]/(np.mean(tmp_eigval)**3)
                    Ani_im_FlA[iix,iiy,iiz]=tmp_eigval[2]/tmp_eigval[1]
                else:
                    Ani_im_RA[iix,iiy,iiz]=np.nan
                    Ani_im_FA[iix,iiy,iiz]=np.nan
                    Ani_im_VA[iix,iiy,iiz]=np.nan
                    Ani_im_FlA[iix,iiy,iiz]=np.nan
                    
                

            dict={"Size box": sbb, "Center box": cbox, "Center mass": center_mass, "Mass": mass, "Eigen value": eigval, "Eigen vector": eigvec,"Main vector": sv3d.setvector3d(mainvec)}
                  
            dict['RA map']=image3d(Ani_im_RA,self.res*size_box)
            dict['FA map']=image3d(Ani_im_FA,self.res*size_box)
            dict['VA map']=image3d(Ani_im_VA,self.res*size_box)
            dict['FlA map']=image3d(Ani_im_FlA,self.res*size_box)
            
            dict['Relative anisotropy'] = np.std(eigval,axis=1)/np.mean(eigval,axis=1)
            dict['Fractional anisotropy'] = np.std(eigval,axis=1)/np.mean(np.multiply(eigval,eigval),axis=1)**0.5
            dict['Volume anisotropy'] = 1.-(np.array(eigval)[:,0]*np.array(eigval)[:,1]*np.array(eigval)[:,2])/(np.mean(eigval,axis=1)**3)
            sortval=np.sort(eigval)
            dict['Flatness anisotropy'] = np.array(sortval)[:,0]/np.array(sortval)[:,1]
            
            list_dict.append(dict)
        
        
        return list_dict
    
    def plotone(self,ax,alpha=0.3):
        id=np.where(self.im==1)
        c=np.ones(len(id[0]))
        ax.scatter(id[0], id[1], id[2], c=c, cmap='viridis', linewidth=0.5,alpha=alpha);
        return ax

                                                 
                                                 
# Function                                               
def tukeywin3D(ss,rad_Tukey):

	iid=np.concatenate((np.linspace(0,rad_Tukey-1,rad_Tukey),np.linspace(ss[0]-rad_Tukey,ss[0]-1,rad_Tukey)))
	jid=np.concatenate((np.linspace(0,rad_Tukey-1,rad_Tukey),np.linspace(ss[1]-rad_Tukey,ss[1]-1,rad_Tukey)))
	kid=np.concatenate((np.linspace(0,rad_Tukey-1,rad_Tukey),np.linspace(ss[2]-rad_Tukey,ss[2]-1,rad_Tukey)))

	rid=np.concatenate((np.linspace(0,rad_Tukey-1,rad_Tukey),np.linspace(rad_Tukey-1,0,rad_Tukey)))

	coeff=np.ones(ss)
	r2=np.zeros(ss)*np.nan
	
	for i in list(range(2*rad_Tukey)):
		r2[int(iid[i]),:,:]=rid[i]
	for i in list(range(2*rad_Tukey)):
		r2[int(rid[i]):ss[0]-int(rid[i]),int(jid[i]),:]=rid[i]
	for i in list(range(2*rad_Tukey)):
		r2[int(rid[i]):ss[0]-int(rid[i]),int(rid[i]):ss[1]-int(rid[i]),int(kid[i])]=+rid[i]

	id=~np.isnan(r2)

	coeff[id]=0.5*(1.+np.cos(np.pi/rad_Tukey*(r2[id]-rad_Tukey)))	

	return coeff

def plotell(center,eigval,eigvec,KT,ax,color='r',alpha=0.7):
    # find the rotation matrix and radii of the axes
    rotation = np.transpose(eigvec)
    radii = (2*KT)**0.5/((eigval)**0.5)

    # calculate cartesian coordinates for the ellipsoid surface
    #print(radii)
    
    u = np.linspace(0.0, 2.0 * np.pi, 60)
    v = np.linspace(0.0, np.pi, 60)
    x = radii[0] * np.outer(np.cos(u), np.sin(v))
    y = radii[1] * np.outer(np.sin(u), np.sin(v))
    z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

    for i in range(len(x)):
        for j in range(len(x)):
            [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center


    ax.plot_surface(x, y, z,  rstride=3, cstride=3,  color=color, linewidth=0.1, alpha=alpha, shade=True)
    
    return

