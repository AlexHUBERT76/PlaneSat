import numpy as np
from astropy.io import ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from PlaneSat import *
import astropy.coordinates as coord

import matplotlib.pyplot as plt




class Galaxy:
    def __init__(self,name):
        self.Name = name
        self.RA = ""#Right ascension astropy icrs format
        self.DE = ""#Declination astropy icrs format
        self.AstroCoord= ""#Skycoordinate with astropy (flatened for satellites)
        self.TrueAstroCoord = ""  # Skycoordinate with astropy
        self.CartCoord = np.array([1,3]) #Cartesian coordinates (x,y,z) flatened at host distance
        self.TrueCartCoord = np.array([1,3]) #Cartesian coordinates (x,y,z)
        self.Type = "" #Giant or dwarf
        self.Dis = "" #Radial distance to Earth in Mpc
        self.Host = ""#Host galaxy
        self.SatNumber = 0
        self.Sat = [] #List of satellites


def SkyPlot(Gal):
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    ax.grid(True)


    ra = []
    dec = []
    Hostra = []
    Hostdec = []
    for i in range(0,len(Gal)-1):

        ra.append(0)
        dec.append(0)
        Hostra.append(0)
        Hostdec.append(0)

        if (Gal[i].Type == "Giant") and (Gal[i].SatNumber != 0):


            #Get galaxy coord
            ra[i] = Gal[i].AstroCoord.ra
            ra[i] = ra[i].wrap_at(180 * u.degree)
            dec[i] = Gal[i].AstroCoord.dec

            ax.scatter(ra[i].radian, dec[i].radian, label = Gal[i].Name)#plotting Giant position
            plt.text(ra[i].radian, dec[i].radian,Gal[i].Name)#plotting Giant name


        elif (Gal[i].Type == "Dwarf") and (Gal[i].Host != ""):

            #Get galaxy coord
            ra[i] = Gal[i].AstroCoord.ra
            ra[i] = ra[i].wrap_at(180 * u.degree)
            dec[i] = Gal[i].AstroCoord.dec

            #Get host galaxy coord
            Hostra[i] = Gal[i].Host.AstroCoord.ra
            Hostra[i] = Hostra[i].wrap_at(180 * u.degree)
            Hostdec[i] = Gal[i].Host.AstroCoord.dec

            plt.plot([ra[i].radian,Hostra[i].radian], [dec[i].radian, Hostdec[i].radian])#Plotting line between satelite and host




    plt.show()

def SatCompil(Gal,hostname,withhost = 0): #Creating catalogue of satelites + host galaxies given atropy catalog Gal and hostname str

    Sat = []
    count = 0
    name = str(hostname)
    for i in range(0,len(Gal)-1):

        if withhost != 0:
            if (Gal[i].Name == name):
                Sat.append(Gal[i])

        if (Gal[i].Host != ""):
            if(Gal[i].Host.Name == name):
                Sat.append(Gal[i])
                count = count + 1

    if Sat == []:
        print("Warning : Data do not contain such host")
    else:
        print(name," has ",count," satellites")

    return Sat

def CartSatCompil(Gal,hostname,dim = 3):

    name = str(hostname)

    Sat = SatCompil(Gal, name)

    Cart = np.zeros([len(Sat)-1, 3])


    if dim == 3:

        for i in range(0, len(Sat) - 1):  # Creating catalogue of satellites cartesian coordinates
            print(Sat[i].Dis)
            Cart[i] = Sat[i].TrueCartCoord

    elif dim == 2:

        for i in range(0, len(Sat) - 1):  # Creating catalogue of satellites cartesian coordinates in a plan
            Cart[i] = Sat[i].CartCoord

    else:

        print("Error : Dimension must be 2 or 3")


    return(Cart)




def main():




    tbl2 = ascii.read("lvg_table2.dat")
    tbl1 = ascii.read("lvg_table1.dat")



    MagTreshold = -20  # minimum absolute magnitude for a galaxy to be an host
    MaxDisHost = 300 #Maximum distance between a dwarf satellite and his host in kpc

    Gal = []
    i = 0



    for row1 in tbl2:



        Gal.append(0)
        Gal[i] = Galaxy(row1["Name"])



        #Formating skycoordinates
        Gal[i].RA = str(row1["RAh"])+"h"+str(row1["RAm"])+"m"+str(row1["RAs"])+"s"
        Gal[i].DE = str(row1["DE-"])+str(row1["DEd"])+"d"+str(row1["DEm"])+"m"+str(row1["DEs"])+"s"
        Gal[i].AstroCoord = SkyCoord(Gal[i].RA, Gal[i].DE, frame='icrs')








        if row1["BMag"] < MagTreshold:  # testing if the galaxy is an host or a dwarf
            Gal[i].Type = "Giant"
        else:
            Gal[i].Type = "Dwarf"

        for row2 in tbl1:


            if row1["Name"] == row2["Name"]: #finding matching name to assign distance for hosts
                Gal[i].Dis = row2["Dis"]
                Gal[i].TrueAstroCoord = SkyCoord(Gal[i].RA, Gal[i].DE, distance=Gal[i].Dis * 10E6 * u.pc, frame='icrs')
                Gal[i].TrueCartCoord = ([Gal[i].TrueAstroCoord.cartesian.x.value, Gal[i].TrueAstroCoord.cartesian.y.value, Gal[i].TrueAstroCoord.cartesian.z.value])
        i = i+1

    SatNumberTotal = 0
    ConflictNumber = 0

    for i in range(0,len(tbl2)-1):#Assigning satellites to host

        if (Gal[i].Type == "Giant") and (Gal[i].Name != "Milky Way"):

            for j in range(0, len(tbl2) - 1):

                sep = Gal[i].TrueAstroCoord.separation_3d(Gal[j].TrueAstroCoord) #Computing separation between host and potential sattelite
                MaxSeparation = MaxDisHost*10E3 *u.pc

                if (sep < MaxSeparation) and (Gal[j].Type == "Dwarf") and (Gal[j].Dis > 0.3 ):

                    #Resolving conflict between two hosts

                    if (Gal[j].Host != "") and (Gal[j].Host.Dis == ""):

                        ConflictNumber = ConflictNumber + 1

                    elif (Gal[j].Host != "") and (Gal[j].Host.Dis != "") and (Gal[j].Host.Dis-Gal[j].Dis < Gal[i].Dis-Gal[j].Dis):

                        dummy = True

                    else:



                        Gal[j].Host = Gal[i]

                        #Adding depth to coordinates
                        Gal[j].AstroCoord = SkyCoord(Gal[j].RA, Gal[j].DE, distance = Gal[j].Host.Dis*10E6*u.pc, frame='icrs')
                        Gal[i].AstroCoord = SkyCoord(Gal[i].RA, Gal[i].DE, distance = Gal[i].Dis*10E6*u.pc, frame='icrs')

                        #Creating list of satellitee for host
                        Gal[i].Sat.append(Gal[j])
                        Gal[i].SatNumber = Gal[i].SatNumber + 1

                        #Assigning cartesian astropy cordinates
                        Gal[j].CartCoord = ([Gal[j].AstroCoord.cartesian.x.value,Gal[j].AstroCoord.cartesian.y.value,Gal[j].AstroCoord.cartesian.z.value])
                        Gal[i].CartCoord = ([Gal[i].AstroCoord.cartesian.x.value,Gal[i].AstroCoord.cartesian.y.value,Gal[i].AstroCoord.cartesian.z.value])







                    SatNumberTotal = SatNumberTotal+1

    HosGal = []
    for i in range(0,len(Gal)-1): #Creating catalogue of host galaxies
        if Gal[i].SatNumber >= 3:
            HosGal.append(Gal[i])
            print(Gal[i].Name,"has", Gal[i].SatNumber, "satellites")

    for i in range(0,len(Gal)-1):
        if Gal[i].SatNumber >= 3:
            HosGal.append(Gal[i])
            print(Gal[i].Name,"has", Gal[i].SatNumber, "satellites")

    CartSatCompil(Gal,"MESSIER081",3)
    SkyPlot(Gal)





if __name__ == '__main__':
    main()