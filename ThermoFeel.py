
"""
ThermoFeel Library
calculates different heat indexes from inputs
"""
#import statements 
import numpy as np


class HeatIndicies():
    """
    Heat Indices is a class that contains all the methods to calculate different heat indexes.
    it includes Relative Humidity, Universal Thermal Climate Index and Mean Radiant Temperature.
    """
    #convert kelvin to celcius
    def kelvin_to_celcius(t):
        t = t - 273.15
        return(t)
    #convert celcius to kelvin
    def celcius_to_kelvin(t):
        t = t + 273.15
        return(t)
    #convert from pa to hpa for e (relative humidity)
    def pa_to_hpa(rh):
        rh = rh/10
        return(rh)
    #Relative Humidity
    def rh(t):
        g = [-2.8365744e3, -6.028076559e3, 1.954263612e1, -2.737830188e-2,
            1.6261698e-5, 7.0229056e-10, -1.8680009e-13, 2.7150305]                      
        tk = t + 273.15
        ess = g[7] * np.log(tk)
        for i in range(7):
            ess = ess + g[i] * pow(tk, (i - 2))  
        ess = np.exp(ess) * 0.01    
        return(ess)
    #Heat Index
    def heat_index(t,rh=None):
        if rh is None:
            rh = HeatIndicies.rh(t)
        t = HeatIndicies.kelvin_to_celcius(t)
        rh = HeatIndicies.pa_to_hpa(rh)
        hiarray = [8.784695, 1.61139411, 2.338549, 0.14611605, 1.2308094 * 10 ** -2, 2.211732 * 10 ** -3, 7.2546 * 10 ** -4, 3.58 * 10 ** -6]
        hi = -hiarray[0] + hiarray[1] * t + hiarray[2] * rh-hiarray[3] * t *\
            rh - hiarray[4] * rh ** 2 + hiarray[5] * t ** 2 * rh+hiarray[6] *\
            t * rh ** 2 - hiarray[7] * t ** 2 * rh ** 2
        return(hi)

    def mean_radiant_temperature(t,ssrd,ssr,fdir,strd,strr,cossza,fp):
        """   
         mrt - Mean Radiant Temperature
         :param T: is 2m temperature
         :param ssrd: is surface solar radiation downwards
         :param ssr: is surface net solar radiation
         :param fdir: is Total sky direct solar radiation at surface
         :param strr: is Surface net thermal radiation
         :param cossza: is Cosine of solar zenith angle
        """
        dsw = ssrd - fdir
        rsw = ssrd - ssr
        lur = strd - strr
        if (cossza > 0.01):
            fdir = fdir / cossza
        mrtcal = 17636684.3 * (0.5 * strd + 0.5 * lur + 0.721649485) \
            * (0.5 * dsw + 0.5 * rsw + fp * fdir)
        mrt = HeatIndicies.kelvin_to_celcius(pow(mrtcal, 0.25))
        return(mrt)
    
    def utci(T,va,mrt=None,rh=None):
        """
        UTCI
        :param t: is 2m temperature
        :param va: is wind speed at 10 meters
        :param mrt: is mean radiant temperature
        :param rh: is relative humidity
        Calculate UTCI with a 6th order polynomial approximation according to:
        Brode, P. et al. Deriving the operational procedure for the 
        Universal Thermal Climate Index (UTCI). Int J Biometeorol (2012) 56: 48.1

        """
        t = HeatIndicies.kelvin_to_celcius(t)
        if rh is None:
            rh = HeatIndicies.rh(t)
        if mrt is None:
            mrt = HeatIndicies.mean_radiant_temperature(t, ssrd, ssr, fdir, strd, strr, cossza, fp)
        e_mrt = np.subtract(mrt, t)
        rh = rh/10.0
          
        if (50 >= t.any() >= -50 and 17 >= va.any() >= 0 and rh.any() <= 5 and 70 >= e_mrt.any() >= -30): 
            t4 = t3 * t
            t3 = t2 * t
            t2 = t * t
            va4 = va3 * va
            va3 = va2 * va
            va2 = va * va
            e_mrt4 = e_mrt3 * e_mrt
            e_mrt3 = e_mrt2 * e_mrt
            e_mrt2 = e_mrt * e_mrt
                
            UTCI_approx = t + 6.07562052E-01 + \
                (-2.27712343E-02) * T +  \
                (8.06470249E-04) * T*T +  \
                (-1.54271372E-04) * T*T*T + \
                (-3.24651735E-06) * T*T*T*T +\
                (7.32602852E-08) * T*T*T*T*T +\
                (1.35959073E-09) * T*T*T*T*T*T +\
                (-2.25836520E+00) * va +  \
                (8.80326035E-02) * T*va +  \
                (2.16844454E-03) * T*T*va + \
                (-1.53347087E-05) * T*T*T*va +\
                (-5.72983704E-07) * T*T*T*T*va +\
                (-2.55090145E-09) * T*T*T*T*T*va +\
                (-7.51269505E-01) * va*va +  \
                (-4.08350271E-03) * T*va*va + \
                (-5.21670675E-05) * T*T*va*va +\
                (1.94544667E-06) * T*T*T*va*va + \
                (1.14099531E-08) * T*T*T*T*va*va +\
                (1.58137256E-01) * va*va*va +  \
                (-6.57263143E-05) * T*va*va*va + \
                (2.22697524E-07) * T*T*va*va*va + \
                (-4.16117031E-08) * T*T*T*va*va*va +\
                (-1.27762753E-02) * va*va*va*va +\
                (9.66891875E-06) * T*va*va*va*va +  \
                (2.52785852E-09) * T*T*va*va*va*va + \
                (4.56306672E-04) * va*va*va*va*va +  \
                (-1.74202546E-07) * T*va*va*va*va*va + \
                (-5.91491269E-06) * va*va*va*va*va*va + \
                (3.98374029E-01) * e_mrt +  \
                (1.83945314E-04) * T*e_mrt + \
                (-1.73754510E-04) * T*T*e_mrt + \
                (-7.60781159E-07) * T*T*T*e_mrt + \
                (3.77830287E-08) * T*T*T*T*e_mrt + \
                (5.43079673E-10) * T*T*T*T*T*e_mrt + \
                (-2.00518269E-02) * va*e_mrt +  \
                (8.92859837E-04) * T*va*e_mrt +  \
                (3.45433048E-06) * T*T*va*e_mrt + \
                (-3.77925774E-07) * T*T*T*va*e_mrt +\
                (-1.69699377E-09) * T*T*T*T*va*e_mrt +\
                (1.69992415E-04) * va*va*e_mrt +  \
                (-4.99204314E-05) * T*va*va*e_mrt + \
                (2.47417178E-07) * T*T*va*va*e_mrt + \
                (1.07596466E-08) * T*T*T*va*va*e_mrt +\
                (8.49242932E-05) * va*va*va*e_mrt +  \
                (1.35191328E-06) * T*va*va*va*e_mrt + \
                (-6.21531254E-09) * T*T*va*va*va*e_mrt +\
                (-4.99410301E-06) * va*va*va*va*e_mrt +  \
                (-1.89489258E-08) * T*va*va*va*va*e_mrt + \
                (8.15300114E-08) * va*va*va*va*va*e_mrt +  \
                (7.55043090E-04) * e_mrt*e_mrt +  \
                (-5.65095215E-05) * T*e_mrt*e_mrt + \
                (-4.52166564E-07) * T*T*e_mrt*e_mrt +\
                (2.46688878E-08) * T*T*T*e_mrt*e_mrt +\
                (2.42674348E-10) * T*T*T*T*e_mrt*e_mrt +\
                (1.54547250E-04) * va*e_mrt*e_mrt +  \
                (5.24110970E-06) * T*va*e_mrt*e_mrt + \
                (-8.75874982E-08) * T*T*va*e_mrt*e_mrt +\
                (-1.50743064E-09) * T*T*T*va*e_mrt*e_mrt +\
                (-1.56236307E-05) * va*va*e_mrt*e_mrt +  \
                (-1.33895614E-07) * T*va*va*e_mrt*e_mrt +  \
                (2.49709824E-09) * T*T*va*va*e_mrt*e_mrt +  \
                (6.51711721E-07) * va*va*va*e_mrt*e_mrt +  \
                (1.94960053E-09) * T*va*va*va*e_mrt*e_mrt + \
                (-1.00361113E-08) * va*va*va*va*e_mrt*e_mrt +\
                (-1.21206673E-05) * e_mrt*e_mrt*e_mrt +  \
                (-2.18203660E-07) * T*e_mrt*e_mrt*e_mrt + \
                (7.51269482E-09) * T*T*e_mrt*e_mrt*e_mrt + \
                (9.79063848E-11) * T*T*T*e_mrt*e_mrt*e_mrt +\
                (1.25006734E-06) * va*e_mrt*e_mrt*e_mrt +  \
                (-1.81584736E-09) * T*va*e_mrt*e_mrt*e_mrt +\
                (-3.52197671E-10) * T*T*va*e_mrt*e_mrt*e_mrt +\
                (-3.36514630E-08) * va*va*e_mrt*e_mrt*e_mrt +  \
                (1.35908359E-10) * T*va*va*e_mrt*e_mrt*e_mrt +  \
                (4.17032620E-10) * va*va*va*e_mrt*e_mrt*e_mrt +  \
                (-1.30369025E-09) * e_mrt*e_mrt*e_mrt*e_mrt +  \
                (4.13908461E-10) * T*e_mrt*e_mrt*e_mrt*e_mrt +  \
                (9.22652254E-12) * T*T*e_mrt*e_mrt*e_mrt*e_mrt + \
                (-5.08220384E-09) * va*e_mrt*e_mrt*e_mrt*e_mrt +  \
                (-2.24730961E-11) * T*va*e_mrt*e_mrt*e_mrt*e_mrt + \
                (1.17139133E-10) * va*va*e_mrt*e_mrt*e_mrt*e_mrt +  \
                (6.62154879E-10) * e_mrt*e_mrt*e_mrt*e_mrt*e_mrt +  \
                (4.03863260E-13) * T*e_mrt*e_mrt*e_mrt*e_mrt*e_mrt + \
                (1.95087203E-12) * va*e_mrt*e_mrt*e_mrt*e_mrt*e_mrt + \
                (-4.73602469E-12) * e_mrt*e_mrt*e_mrt*e_mrt*e_mrt*e_mrt +\
                (5.12733497E+00) * rh +  \
                (-3.12788561E-01) * T*rh + \
                (-1.96701861E-02) * T*T*rh + \
                (9.99690870E-04) * T*T*T*rh + \
                (9.51738512E-06) * T*T*T*T*rh + \
                (-4.66426341E-07) * T*T*T*T*T*rh +\
                (5.48050612E-01) * va*rh +  \
                (-3.30552823E-03) * T*va*rh + \
                (-1.64119440E-03) * T*T*va*rh + \
                (-5.16670694E-06) * T*T*T*va*rh + \
                (9.52692432E-07) * T*T*T*T*va*rh + \
                (-4.29223622E-02) * va*va*rh +  \
                (5.00845667E-03) * T*va*va*rh +  \
                (1.00601257E-06) * T*T*va*va*rh + \
                (-1.81748644E-06) * T*T*T*va*va*rh +\
                (-1.25813502E-03) * va*va*va*rh +  \
                (-1.79330391E-04) * T*va*va*va*rh + \
                (2.34994441E-06) * T*T*va*va*va*rh + \
                (1.29735808E-04) * va*va*va*va*rh +  \
                (1.29064870E-06) * T*va*va*va*va*rh + \
                (-2.28558686E-06) * va*va*va*va*va*rh + \
                (-3.69476348E-02) * e_mrt*rh +  \
                (1.62325322E-03) * T*e_mrt*rh +  \
                (-3.14279680E-05) * T*T*e_mrt*rh + \
                (2.59835559E-06) * T*T*T*e_mrt*rh + \
                (-4.77136523E-08) * T*T*T*T*e_mrt*rh +\
                (8.64203390E-03) * va*e_mrt*rh +  \
                (-6.87405181E-04) * T*va*e_mrt*rh +\
                (-9.13863872E-06) * T*T*va*e_mrt*rh +\
                (5.15916806E-07) * T*T*T*va*e_mrt*rh + \
                (-3.59217476E-05) * va*va*e_mrt*rh +  \
                (3.28696511E-05) * T*va*va*e_mrt*rh +  \
                (-7.10542454E-07) * T*T*va*va*e_mrt*rh +\
                (-1.24382300E-05) * va*va*va*e_mrt*rh +  \
                (-7.38584400E-09) * T*va*va*va*e_mrt*rh + \
                (2.20609296E-07) * va*va*va*va*e_mrt*rh +  \
                (-7.32469180E-04) * e_mrt*e_mrt*rh +  \
                (-1.87381964E-05) * T*e_mrt*e_mrt*rh + \
                (4.80925239E-06) * T*T*e_mrt*e_mrt*rh + \
                (-8.75492040E-08) * T*T*T*e_mrt*e_mrt*rh +\
                (2.77862930E-05) * va*e_mrt*e_mrt*rh +  \
                (-5.06004592E-06) * T*va*e_mrt*e_mrt*rh +\
                (1.14325367E-07) * T*T*va*e_mrt*e_mrt*rh + \
                (2.53016723E-06) * va*va*e_mrt*e_mrt*rh +  \
                (-1.72857035E-08) * T*va*va*e_mrt*e_mrt*rh +\
                (-3.95079398E-08) * va*va*va*e_mrt*e_mrt*rh +\
                (-3.59413173E-07) * e_mrt*e_mrt*e_mrt*rh +  \
                (7.04388046E-07) * T*e_mrt*e_mrt*e_mrt*rh +  \
                (-1.89309167E-08) * T*T*e_mrt*e_mrt*e_mrt*rh +\
                (-4.79768731E-07) * va*e_mrt*e_mrt*e_mrt*rh +  \
                (7.96079978E-09) * T*va*e_mrt*e_mrt*e_mrt*rh +  \
                (1.62897058E-09) * va*va*e_mrt*e_mrt*e_mrt*rh +  \
                (3.94367674E-08) * e_mrt*e_mrt*e_mrt*e_mrt*rh +  \
                (-1.18566247E-09) * T*e_mrt*e_mrt*e_mrt*e_mrt*rh +\
                (3.34678041E-10) * va*e_mrt*e_mrt*e_mrt*e_mrt*rh + \
                (-1.15606447E-10) * e_mrt*e_mrt*e_mrt*e_mrt*e_mrt*rh +\
                (-2.80626406E+00) * rh*rh +  \
                (5.48712484E-01) * T*rh*rh +  \
                (-3.99428410E-03) * T*T*rh*rh +\
                (-9.54009191E-04) * T*T*T*rh*rh +\
                (1.93090978E-05) * T*T*T*T*rh*rh + \
                (-3.08806365E-01) * va*rh*rh +  \
                (1.16952364E-02) * T*va*rh*rh +  \
                (4.95271903E-04) * T*T*va*rh*rh + \
                (-1.90710882E-05) * T*T*T*va*rh*rh +\
                (2.10787756E-03) * va*va*rh*rh +  \
                (-6.98445738E-04) * T*va*va*rh*rh +\
                (2.30109073E-05) * T*T*va*va*rh*rh +\
                (4.17856590E-04) * va*va*va*rh*rh +  \
                (-1.27043871E-05) * T*va*va*va*rh*rh +\
                (-3.04620472E-06) * va*va*va*va*rh*rh +\
                (5.14507424E-02) * e_mrt*rh*rh +  \
                (-4.32510997E-03) * T*e_mrt*rh*rh +\
                (8.99281156E-05) * T*T*e_mrt*rh*rh +\
                (-7.14663943E-07) * T*T*T*e_mrt*rh*rh +\
                (-2.66016305E-04) * va*e_mrt*rh*rh +  \
                (2.63789586E-04) * T*va*e_mrt*rh*rh +  \
                (-7.01199003E-06) * T*T*va*e_mrt*rh*rh + \
                (-1.06823306E-04) * va*va*e_mrt*rh*rh +  \
                (3.61341136E-06) * T*va*va*e_mrt*rh*rh +  \
                (2.29748967E-07) * va*va*va*e_mrt*rh*rh +  \
                (3.04788893E-04) * e_mrt*e_mrt*rh*rh +  \
                (-6.42070836E-05) * T*e_mrt*e_mrt*rh*rh +\
                (1.16257971E-06) * T*T*e_mrt*e_mrt*rh*rh + \
                (7.68023384E-06) * va*e_mrt*e_mrt*rh*rh +  \
                (-5.47446896E-07) * T*va*e_mrt*e_mrt*rh*rh +\
                (-3.59937910E-08) * va*va*e_mrt*e_mrt*rh*rh +\
                (-4.36497725E-06) * e_mrt*e_mrt*e_mrt*rh*rh + \
                (1.68737969E-07) * T*e_mrt*e_mrt*e_mrt*rh*rh + \
                (2.67489271E-08) * va*e_mrt*e_mrt*e_mrt*rh*rh + \
                (3.23926897E-09) * e_mrt*e_mrt*e_mrt*e_mrt*rh*rh +\
                (-3.53874123E-02) * rh*rh*rh +  \
                (-2.21201190E-01) * T*rh*rh*rh + \
                (1.55126038E-02) * T*T*rh*rh*rh + \
                (-2.63917279E-04) * T*T*T*rh*rh*rh +\
                (4.53433455E-02) * va*rh*rh*rh +  \
                (-4.32943862E-03) * T*va*rh*rh*rh +\
                (1.45389826E-04) * T*T*va*rh*rh*rh +\
                (2.17508610E-04) * va*va*rh*rh*rh +  \
                (-6.66724702E-05) * T*va*va*rh*rh*rh +\
                (3.33217140E-05) * va*va*va*rh*rh*rh + \
                (-2.26921615E-03) * e_mrt*rh*rh*rh +  \
                (3.80261982E-04) * T*e_mrt*rh*rh*rh +  \
                (-5.45314314E-09) * T*T*e_mrt*rh*rh*rh +\
                (-7.96355448E-04) * va*e_mrt*rh*rh*rh +  \
                (2.53458034E-05) * T*va*e_mrt*rh*rh*rh +  \
                (-6.31223658E-06) * va*va*e_mrt*rh*rh*rh + \
                (3.02122035E-04) * e_mrt*e_mrt*rh*rh*rh +  \
                (-4.77403547E-06) * T*e_mrt*e_mrt*rh*rh*rh +\
                (1.73825715E-06) * va*e_mrt*e_mrt*rh*rh*rh + \
                (-4.09087898E-07) * e_mrt*e_mrt*e_mrt*rh*rh*rh +\
                (6.14155345E-01) * rh*rh*rh*rh +  \
                (-6.16755931E-02) * T*rh*rh*rh*rh +\
                (1.33374846E-03) * T*T*rh*rh*rh*rh +\
                (3.55375387E-03) * va*rh*rh*rh*rh +  \
                (-5.13027851E-04) * T*va*rh*rh*rh*rh +\
                (1.02449757E-04) * va*va*rh*rh*rh*rh + \
                (-1.48526421E-03) * e_mrt*rh*rh*rh*rh + \
                (-4.11469183E-05) * T*e_mrt*rh*rh*rh*rh +\
                (-6.80434415E-06) * va*e_mrt*rh*rh*rh*rh +\
                (-9.77675906E-06) * e_mrt*e_mrt*rh*rh*rh*rh +\
                (8.82773108E-02) * rh*rh*rh*rh*rh +  \
                (-3.01859306E-03) * T*rh*rh*rh*rh*rh +\
                (1.04452989E-03) * va*rh*rh*rh*rh*rh + \
                (2.47090539E-04) * e_mrt*rh*rh*rh*rh*rh +\
                (1.48348065E-03) * rh*rh*rh*rh*rh*rh
            return(UTCI_approx)
                
        
        
        
#The main method to run the code
if __name__ == '__main__':
    print("main") 
     #my_example_nc_file2 = 'C:/Users/moc2/OneDrive - University of Reading/projaveragetemp.nc'
     # fh2 = Dataset(my_example_nc_file2, mode='r')
      #lons = fh2.variables['longitude'][:]
      #lats = fh2.variables['latitude'][:]
      #T=fh2.variables['layer'][:]
     # heatindex = HeatIndicies.heat_index(T)
     # utcitest = HeatIndicies.utci(T, va=T, mrt=T, rh=T)
     # print(utcitest)
     # print(heatindex)
     # lon, lat = np.meshgrid(lons, lats)
     # fig = plt.figure()
     # ax = plt.axes(projection=ccrs.PlateCarree())
     # ax.coastlines()
     # filled_c=plt.contourf(lons, lats, heatindex, 60,transform=ccrs.PlateCarree(),cmap="YlOrRd",clim=(-50.50))
     # fig.colorbar(filled_c, orientation='horizontal',ticks=[-50,-40,-30,-20,-10,0,10,20,30,40,50])
     # plt.show()


     