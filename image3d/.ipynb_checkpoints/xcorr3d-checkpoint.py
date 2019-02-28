import image3d.image3d as im3d
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import image3d.uniform_dist as uniform_dist
import setvector3d.setvector3d as sv3d
import matplotlib.tri as tri
import skimage.feature as skf

class xcorr3d(im3d.image3d):
    '''
    3D Autocorrelation function
    '''

    pass

    def __init__(self,data,res,Cinf):
        '''
		Construtor of xcorr3d object
		:param data: array Nx3 dimention
		:type data: np.array
		:param res: resolution of the image (default 1)
		:type res: float
		:param Cinf:
		:type Cinf: float
		'''
        
        self.im=data
        self.res=res
        self.Cinf=Cinf
        return
    
    def local_maxima(self,nb=20,win=20):
        '''
        Detect the local maximum of the autocorrelation function
        :param th: threshold percentage in the range [Cinf 1]
        :type th: float
        :param win: windows size for filter NxNxN
        :type win: int
        '''
                
        coordinates = skf.peak_local_max(self.im, min_distance=win,indices=True,exclude_border=False,num_peaks=nb,threshold_abs=self.Cinf)
        
        ss=np.shape(self.im)
        vec=np.zeros(np.shape(coordinates))
        vec[:,0]=coordinates[:,0]-ss[0]/2
        vec[:,1]=coordinates[:,1]-ss[1]/2
        vec[:,2]=coordinates[:,2]-ss[2]/2
        id=np.where((vec[:,2]>0)*(np.sum(~np.isnan(vec),axis=1)>0) )
        
        vc=sv3d.setvector3d(vec[id[0],:])
        coordinates=coordinates[id[0],:]
        val=np.zeros(coordinates.shape[0])
        for i in list(range(coordinates.shape[0])):
            val[i]=self.im[coordinates[i,0],coordinates[i,1],coordinates[i,2]]
        
        return coordinates,vc,val
        
    
    
    def PhiTheta_corr_length(self):
        '''
		Extract the correlation length in function of Theta and phi
		Theta is the angle in the plan xOy between [0,2*Pi]
		phi is the angle from z between [0,Pi/2]
		'''
        Theta=np.arange(360)*np.pi/180.
        Phi=np.arange(90)*np.pi/180.
        ss=np.shape(self.im)

		#center=int(ss[0]/2.)
        c0=int(ss[0]/2.)
        c1=int(ss[1]/2.)
        c2=int(ss[2]/2.)

		#rmax=np.max(ss)/2.-1 #here the minus 1 is only to not overshoot into the matrix	
        a=np.zeros(3)
        a[0]=ss[0]/2.-1.
        a[1]=ss[1]/2.-1.
        a[2]=ss[2]/2.-1.	

        #init=np.ones(3)*center
        init=np.array([c0,c1,c2])
		
        xmin=np.zeros([360,90])
        ki=0
        for i in Theta:
            kj=0
            for j in Phi:
                b=np.array([abs(np.sin(j)*np.cos(i)),abs(np.sin(j)*np.sin(i)),abs(np.cos(j))])
                idb=np.where(b!=0)
                rmax=np.min(a[idb[0]]/b[idb[0]])
                
                xe=init[0]+np.sin(j)*np.cos(i)*rmax
                ye=init[1]+np.sin(j)*np.sin(i)*rmax
                ze=init[2]+np.cos(j)*rmax
                
                end=np.array([xe,ye,ze])
				#extract data on the line Theta Phi
                [res,xl]=self.extract_profil(init,end)
				#find where the value of Autocorrelation function is inferior to Cinf
                id=np.where(res<self.Cinf)
                if np.size(id)==0:
                    xmin[ki,kj]=np.inf
                else:
                    xmin[ki,kj]=xl[id[0][0]]
                kj=kj+1
            ki=ki+1
            
        plt.imshow(xmin)
        plt.axis('equal')
        plt.colorbar(orientation='vertical',aspect=4,shrink=0.5)
        return
    
    def stereographic_corr_length(self,output='No',coeffCinf=np.array([1]),pc=10,usePI=False):
        '''
        :param output: Destination of the figure file
        :type output: str
        :param coeffCinf: array of coeeficient you want to applied to Cinf (dimention N, means N figure plotted)
        :type coeffCinf: np.array
        :param pc: percentage of highest orientation taken in the statistic (Between 0 and 100)
        :type pc: float
        :param usePI: compute the autocorrelation radius using the "Porte Integrale" (default False) overwise It compute the intersection between Cinf and Ax
        :type usePI: bool
        '''
        val=uniform_dist.unidist
        dim=int(np.size(val)/3)
        xx=val.reshape([dim,3])
        id=np.where(xx[:,2]>0)
        # unit vector where I compute the correlation length
        xxuse=xx[id[0],:]
        vuse=sv3d.setvector3d(xxuse)
        # add point on the disc for contourf
        tot=1000
        omega = np.linspace(0, 2*np.pi, tot)
        zcir = np.zeros(tot)
        xcir = np.cos(omega)
        ycir = np.sin(omega)
        vcir=sv3d.setvector3d(np.transpose(np.array([xcir,ycir,zcir])))
            
        vtot=vuse.concatenate(vcir)
        
        #############################
        ## Compute the pole figure ##
        #############################
        # create initial point
        ss=np.shape(self.im)
        #center=int(ss[0]/2.)
        c0=int(ss[0]/2.)
        c1=int(ss[1]/2.)
        c2=int(ss[2]/2.)
        #rmax=np.max(ss)/2.-1 #here the minus 1 is only to not overshoot into the matrix	
        a=np.zeros(3)
        a[0]=ss[0]/2.-1.
        a[1]=ss[1]/2.-1.
        a[2]=ss[2]/2.-1.	
		
        init=np.array([c0,c1,c2])
        nb_points=np.shape(vtot.vector)[0]
        nbimg=np.size(coeffCinf)
        xmin=np.zeros([nb_points,nbimg])
        rmax=np.zeros(nb_points)
        
		
        for i in list(range(nb_points)): # loop for the differente direction in autocorrelation function
            b=np.array([abs(vtot.vector[i,0]),abs(vtot.vector[i,1]),abs(vtot.vector[i,2])])
            idb=np.where(b!=0)
            rmax[i]=np.min(a[idb[0]]/b[idb[0]])
            xe=init[0]+vtot.vector[i,0]*rmax[i]
            ye=init[1]+vtot.vector[i,1]*rmax[i]
            ze=init[2]+vtot.vector[i,2]*rmax[i]
            end=np.array([xe,ye,ze])
            [res,xl]=self.extract_profil(init,end)
            
            for j in list(range(nbimg)): # loop for the differente value of Cinf
                if usePI:
                    xmin[i,j]=2.*np.trapz(res-self.Cinf*coeffCinf[j],xl)
                else:
                    id=np.where(res < self.Cinf*coeffCinf[j])
                    if np.size(id)==0:
                        xmin[i,j]=np.inf
                    else:
                        xmin[i,j]=xl[id[0][0]]
        #############################           
        ### Compute the statistics ##
        #############################
        # Stat on correlation length #
        radiusmean=np.zeros(nbimg)  # Mean value of the correlation radius for pc% highest point
        radiusstd=np.zeros(nbimg)   # Standard deviation value of the correlation radius for pc% highest point
        # Stat using second order orientation tensor on set of orientation with thes pc% highest value of correlation radius
        eigvector=np.zeros([nbimg,3,3]) # eigvector[j,:,i] correspond to w[j,i] where j is the increment on different Cinf value
        eigvalue=np.zeros([nbimg,3]) #
        
        
        for j in list(range(nbimg)): # loop for the differente value of Cinf
            val=np.nanpercentile(xmin[:,j],100.-pc) # Find the limite for correlation raidus
            if np.isinf(val):
                id=np.where(np.isinf(xmin[0:nb_points-tot,j]))
            else:
                id=np.where(xmin[0:nb_points-tot,j]>val) # Find the position for the orientation higher than the correlation radius. I remove the evaluation on the cicle that I add to be sure that a sampling my sphere hogeneusly
            vmax=sv3d.setvector3d(vtot.vector[id[0],:])
            #vmax.stereoplot()
            # statistic about radius correlation
            radiusmean[j]=np.nanmean(vtot.vector[id[0],j])
            radiusstd[j]=np.nanstd(vtot.vector[id[0],j])
            
            #############################################################################
            # Compute the second order orientation tensor for the subseb of orientation #
            #############################################################################
            # The good way to do it
            #Tensor=np.zeros([3,3])
            #for i in list(range(id[0])):
            #    Tensor=Tensor+1./len(id[0])*np.multiply.outer(xxuse[id[0]],xxuse[id[0]])

            eigvalue[j,:],eigvector[j,:,:]=vmax.OrientationTensor2nd()
            
        #####################
        ## Build the Image ##
        #####################
        # Compute for image
        # Make all the vector pointing in the upper hemisphere
        phip,thetap=vtot.cart2spher()
        xx = np.multiply(2*np.sin(phip/2),np.cos(thetap))
        yy = np.multiply(2*np.sin(phip/2),np.sin(thetap))
        # Value to plot circle
        angle=np.array([90.,60.,30.])
        rci=2.*np.sin(angle/2.*np.pi/180.)
        
             
        for j in list(range(nbimg)):
            plt.figure(figsize=(10,10),dpi=160)
            triang = tri.Triangulation(xx, yy)
            zz=xmin[:,j]
            mm=np.max(zz[~np.isinf(zz)])
            zz[np.isinf(zz)]=mm
            ##### Deal with inf value
            plt.tricontourf(triang, zz, 10)
            plt.colorbar(orientation='vertical',aspect=4,shrink=0.5)
            if len(np.where(zz==mm)[0])>1:
                idd=np.where(zz==mm)
                zzm=np.ones(zz.shape)
                zzm[idd]=0                
                plt.tricontourf(triang, zzm, 10,cmap=cm.binary_r)
                plt.clim(1,1)
                
                
            
            #idinf=np.isinf(xmin[:,j])
            #plt.contourf(xx[~idinf],yy[~idinf], xmin[~idinf,j], 10)
            #plt.contourf(xx[idinf],yy[idinf], -np.ones(np.sum(idinf)), 10)
            ### STOP HERE ####
            #plt.clim(np.min(Zc[:])-1., np.max(Zc[:])+1.)
           
            # compute a 3 circle
            
            
            
            for i in list(range(len(rci))):
                omega = np.linspace(0, 2*np.pi, 1000)
                x_circle = rci[i]*np.cos(omega)
                y_circle = rci[i]*np.cos(i*np.pi/180.)*np.sin(omega)
                if i==0:
                    plt.plot(x_circle, y_circle,'k', linewidth=3)
                else:
                    plt.plot(x_circle, y_circle,'k', linewidth=1.5)
                    plt.text(x_circle[200], y_circle[200]+0.04,'$\phi$='+str(angle[i])+'°')
            
            # plot Theta line
            plt.plot([0,0],[-1*rci[0],1*rci[0]],'k', linewidth=1.5)
            plt.text(rci[0]-0.2, 0+0.06,'$\Theta$=0°')
            plt.text(-rci[0]+0.1, 0-0.06,'$\Theta$=180°')
            plt.plot([-rci[0],rci[0]],[0,0],'k', linewidth=1.5)
            plt.text(-0.25, rci[0]-0.25,'$\Theta$=90°')
            plt.text(0.01, -rci[0]+0.15,'$\Theta$=270°')
            plt.plot([-0.7071*rci[0],0.7071*rci[0]],[-0.7071*rci[0],0.7071*rci[0]],'k', linewidth=1.5)
            plt.plot([-0.7071*rci[0],0.7071*rci[0]],[0.7071*rci[0],-0.7071*rci[0]],'k', linewidth=1.5)
            
            
            # draw a cross for x and y direction
            plt.plot([1*rci[0], 0],[0, 1*rci[0]],'+k',markersize=12)
            # write axis
            plt.text(1.05*rci[0], 0, r'X')
            plt.text(0, 1.05*rci[0], r'Y')
            plt.axis('equal')
            plt.axis('off')
            ########################
            # plot the eigen value #
            ########################
            for i in list(range(3)): # Loop on the 3 eigenvalue
                if (eigvector[j,2,i]<0):
                    v=-eigvector[j,:,i]
                else:
                    v=eigvector[j,:,i]
                    
                phiee=np.arccos(v[2])
                thetaee=np.arctan2(v[1],v[0])
                xxv = np.multiply(2*np.sin(phiee/2),np.cos(thetaee))
                yyv = np.multiply(2*np.sin(phiee/2),np.sin(thetaee))
                    
                plt.plot(xxv,yyv,'sk',markersize=8)
                plt.text(xxv+0.04, yyv+0.04,str(round(eigvalue[j,i],2)))
                    
                        
            if output != 'No':
                plt.xlabel('Stereoprojection - radius length')
                plt.savefig(output + '_Cinfx ' + str(coeffCinf[j]) + '.png')
                
                
                
        return eigvalue,eigvector,radiusmean,radiusstd,xmin,rmax
    
    
    
    
    def correlation_profil(self,vector):
        ''' 
        Plot the correlation profile for a given direction vector
        :param vector: unit vector length in 3 dimention
        :type vector: np.array
        '''
        if np.linalg.norm(vector)!=1:
            vector=vector/np.linalg.norm(vector)
            print('Normalized vector to 1')
            
        ss=np.shape(self.im)
        
        c0=int(ss[0]/2.)
        c1=int(ss[1]/2.)
        c2=int(ss[2]/2.)
        
        a=np.zeros(3)
        a[0]=ss[0]/2.-1.
        a[1]=ss[1]/2.-1.
        a[2]=ss[2]/2.-1.
        # intial point.
        init=np.array([c0,c1,c2])
        
        b=np.array([abs(vector[0]),abs(vector[1]),abs(vector[2])])
        idb=np.where(b!=0)
        
        rmax=np.min(a[idb[0]]/b[idb[0]])
        
        xe=init[0]+vector[0]*rmax
        ye=init[1]+vector[1]*rmax
        ze=init[2]+vector[2]*rmax
        end=np.array([xe,ye,ze])
        [res,xl]=self.extract_profil(init,end)
            
        plt.figure()
        plt.plot(xl,res)

        plt.plot([xl[0],xl[-1]],[self.Cinf,self.Cinf])
        plt.xlabel('Distance')
        plt.ylabel('Auto-correlation value')
        
        return
