import numpy as np
import matplotlib.pyplot as plt 


class MobiusStrip:
    def __init__(self, R, w, n): #constructor
        self.R = R  # radius
        self.w = w  # width
        self.n = n  # mesh density

        # Generate the grid of parameters u and v
        self.u = np.linspace(0, 2*np.pi, self.n) #linspace will makesure that u lies b/w [0,2pi], n evenly spaced values of v.
        self.v = np.linspace(-self.w/2, self.w/2, self.n) # n evenly spaced values of v.
        self.U, self.V = np.meshgrid(self.u, self.v) #returns two 2D arrays consisting of values of grid coordinates u[i,j],v[i,j], U of size nxn and V of size nxn

       
        self.X, self.Y, self.Z = self._generate_mesh() #returns the co-ordinates of MobiusStrip

    def _generate_mesh(self):
        #Computes the 3D coordinates of the Mobius strip.
        U = self.U
        V = self.V
        X = (self.R + V*np.cos(U/2))*np.cos(U) 
        Y = (self.R + V*np.cos(U/2))*np.sin(U)
        Z = V * np.sin(U/2)
        return X, Y, Z

    def surface_area(self):
        """
        Numerically approximates the surface area using the magnitude of the cross product
        of partial derivatives with respect to u and v.
        """
        du = 2*np.pi/(self.n-1)     # both du and dv used to approximate the value of surfave area 
        dv = self.w/(self.n-1)

        Xu = np.gradient(self.X, axis=1) / du    #Approximate partial derivative with respect to du
        Yu = np.gradient(self.Y, axis=1) / du
        Zu = np.gradient(self.Z, axis=1) / du

        Xv = np.gradient(self.X, axis=0) / dv    #Approximate partial derivative with respect to dv
        Yv = np.gradient(self.Y, axis=0) / dv
        Zv = np.gradient(self.Z, axis=0) / dv

        # Cross product components
        Nx = Yu*Zv-Zu*Yv  
        Ny = Zu*Xv-Xu*Zv
        Nz = Xu*Yv-Yu*Xv
        
        dA = np.sqrt(Nx**2 + Ny**2 + Nz**2) #computes the magnitude of the vectors
        return np.sum(dA) * du * dv # returns by summing all small surface areas

    def edge_length(self):
        """
        Approximates the length of the boundary curve.
        The boundary consists of two edges: v = -w/2 and v = w/2.
        The total edge length is the sum of these.
        """
        edges = []
        for edge_v in [-self.w/2, self.w/2]:
            u = self.u
            x = (self.R+edge_v*np.cos(u/2))*np.cos(u)
            y = (self.R+edge_v*np.cos(u/2))*np.sin(u)
            z = edge_v*np.sin(u/2)
            dx = np.diff(x)
            dy = np.diff(y)
            dz = np.diff(z)
            edge_length = np.sum(np.sqrt(dx**2 + dy**2 + dz**2)) 
            edges.append(edge_length)
        return sum(edges)

    def plot(self):
        """Generates a 3D surface plot of the Mobius strip."""
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(self.X, self.Y, self.Z, color='blue', edgecolor='k', linewidth=0.5, rstride=1, cstride=1) 
        ax.set_title("Möbius Strip")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.tight_layout()
        plt.show()
        
        
if __name__ == "__main__":
    R = 5
    w = 1
    n = 300

    mobius = MobiusStrip(R, w, n)
    area = mobius.surface_area()
    edge = mobius.edge_length()

    print(f"Approximated Surface Area: {area:.4f}")
    print(f"Approximated Edge Length: {edge:.4f}")
    mobius.plot() 
