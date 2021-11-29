import numpy as np
import scipy.stats as ss

class RV:
    def __init__(self, dist_name: str,loc: float, scale: float):
        '''
        dist_name: 'norm' or 'uniform'
        loc, scale: the same as used in scipy.stats
        '''
        
        self.dist_name = dist_name
        self.loc = loc
        self.scale = scale
        
        if self.dist_name == 'uniform': 
            self.rvf = np.random.uniform
            self.rvf_scale = self.loc + self.scale
            self.pdf = ss.uniform.pdf
        
            
        if self.dist_name == 'normal':
            self.rvf = np.random.normal
            self.rvf_scale = self.scale
            self.pdf = ss.norm.pdf


    def __lt__(self, other: 'RV'):
        return self.loc < other.loc
        

    def rvfv(self):
        '''
        Generate a random sample from the given prior distribution
        '''
        return self.rvf(self.loc,self.rvf_scale)
    
    def pdfv(self,x):
        '''
        Get the pdf value for a give value x
        '''
        return self.pdf(x,self.loc,self.scale)
    
    def mutate(self):
        '''
        Generate a random sample from the given prior distribution and use it as a the location for a new RV object
        '''
        return RV(dist_name=self.dist_name, loc=self.rvfv(), scale=self.scale)

