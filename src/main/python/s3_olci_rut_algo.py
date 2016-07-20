# -*- coding: utf-8 -*-
"""
Created on Thurs May 26 17:51:33 2016

@author: shunt
"""

import numpy as np
import s3_olci_l1_rad_conf as rad_conf


class S3OLCIRutAlgo:
    """
    Algorithm for the Sentinel-3 OLCI Radiometric Uncertainty Tool (RUT)
    """

    def __init__(self):
        # uncertainty values for DS and Abs.cal
        self.u_diffchar_2 = rad_conf.u_diffchar_2
        self.u_diffmod_2 = rad_conf.u_diffmod_2
        self.u_diffalign_2 = rad_conf.u_diffalign_2
        self.u_calistr_2 = rad_conf.u_calistr_2
        self.u_postcalistr_2 = rad_conf.u_postcalistr_2
        self.u_calispeck_2 = rad_conf.u_calispeck_2
        
        self.k = 1

        self.unc_select = [True, True, True, True, True, True, True, True, True, True, True,
                           True]  # list of booleans with user selected uncertainty sources(order as in interface)

    def unc_calculation(self, band_data, band_id):
        """
        This function represents the core of the S3 OLCI RUTv1.
        It takes as an input the pixel data of a specific band and tile in
        a S3-OLCI-L1B product and produces an image with the same dimensions that
        contains the radiometric uncertainty of each pixel radiance.

        :param band_data: list with the quantized L1B radiance pixels of a band (flattened; 1-d)
        :param band_id: zero-based index of the band
        :return: list of u_int8 with uncertainty associated to each pixel.
        """

        #######################################################################
        # 1.	Initial check
        #######################################################################
        # to do

        #######################################################################
        # 2.	Instrument Simulator
        #######################################################################
        
        #OLCIRUTv2.0 would include way to differentiate modules - for now use centre as representative of whole instrument
        mod_id = 3
        
        # Estimate Level 0 counts from input Level 1
        X = band_data * rad_conf.C_RSs[mod_id, band_id] # propagate to CCD

        # Open dark signal estimate values
        X_DS = rad_conf.DSs[mod_id, band_id]
        
        # Estimate smear band value
        # (V1.0 of S3-OLCI-RUT Planned to make use of better calc_smear_band()
        # not enough time to fully implement in this version)
        X_sm = rad_conf.X_sm_L_ref/rad_conf.X_ref*X

        X_tot = [(X + X_sm + DS)*n_pxl for DS, X, n_pxl in zip(X_DS, X, rad_conf.n_pxls[mod_id,band_id])]/rad_conf.digi_steps[mod_id, band_id]
        X_sm = X_sm/rad_conf.digi_steps_sm[mod_id]

        #######################################################################        
        # 3.	L0 uncertainty contributors: noise
        #######################################################################

        if self.unc_select[0]:
            u_noise_2 = X_tot**2 + rad_conf.a_noise[band_id]**2
        else:
            u_noise_2 = 0

        if self.unc_select[1]:
            u_instage = rad_conf.u_instage[band_id]
        else:
            u_instage = 0    

        if self.unc_select[1]:
            u_ccdstab_2 = rad_conf.u_ccdstab_2[band_id]
        else:
            u_ccdstab_2 = 0    

        if self.unc_select[2]:
            u_PS = (rad_conf.a_PS * X**-2 + rad_conf.b_PS)**0.5
        else:
            u_PS = 0 
        
        #######################################################################
        # 4.    L1B uncertainty contributors: absolute calibration coefficient
        #######################################################################
        
        if not self.unc_select[3]:
            self.u_diffchar_2 = 0  # predefined but updated to 0 if deselected by user

        if not self.unc_select[4]:
            self.u_diffmod_2 = 0  # predefined but updated to 0 if deselected by user
        
        if not self.unc_select[5]:
            self.u_diffalign = 0  # predefined but updated to 0 if deselected by user
        
        if not self.unc_select[6]:
            self.u_calistr_2 = 0  # predefined but updated to 0 if deselected by user
            
        if not self.unc_select[7]:
            self.u_postcalistr_2 = 0  # predefined but updated to 0 if deselected by user
            
        if not self.unc_select[8]:
            self.u_calispeck_2 = 0  # predefined but updated to 0 if deselected by user
            
        if self.unc_select[9]:
            u_diff1age = rad_conf.u_diff1age[band_id]
        else:
            u_diff1age = 0
            
        if self.unc_select[10]:
            u_diff2age = rad_conf.u_diff2age[band_id]
        else:
            u_diff2age = 0
        
        #######################################################################        
        # 5.	L1B uncertainty contributors: non-linarity correction
        #######################################################################        

        if self.unc_select[11]:
            u_INL_2 = rad_conf.a_INL(rad_conf.b_INL[band_id]((X_sm/X)**2+rad_conf.c_INL[band_id])+((X_tot/X)**2+rad_conf.d_INL[band_id]))
        else:
            u_INL_2 = 0
            
        if self.unc_select[12]:
            u_DNL_2 = rad_conf.a_DNL[band_id] * X**-2  + rad_conf.b_DNL[band_id]
        else:
            u_DNL_2 = 0

        #######################################################################
        # 6.	L1B uncertainty contributors: dark signal correction
        #######################################################################

        if self.unc_select[7]:
            u_off_2 = rad_conf.a_off * X**-2 + rad_conf.b_off
        else:
            u_off_2 = 0
            
        if self.unc_select[7]:
            u_darkstab_2 = rad_conf.a_darkstab * X**-2
        else:
            u_darkstab_2 = 0

        #######################################################################
        # 7.    L1B uncertainty contributors: smear correction
        #######################################################################

        if self.unc_select[7]:
            u_SGR_2 = rad_conf.a_SGR((X_sm/X)+rad_conf.b_SGR)**2+rad_conf.c_SGR((X_sm/X)**2+rad_conf.b_SGR**2)
        else:
            u_SGR_2 = 0

        #######################################################################
        # 8.    L1B uncertainty contributors: stray light correction
        #######################################################################

        # TBC

        #######################################################################        
        # 9.	Combine uncertainty contributors
        #######################################################################        
        # values given as percentages. Multiplied by 10 and saved to 1 byte(uint8)
        # Clips values to 0-250 --> uncertainty >=25%  assigns a value 250.
        
        u_1sigma = u_instage + u_diff1age + u_diff2age + u_PS + (self.u_diffchar_2 + \
                   self.u_diffmod_2 + self.u_diffalign_2 + u_ccdstab_2 + \
                   self.u_calistr_2 + self.u_postcalistr_2 + self.u_calispec_2 + u_off_2 + \
                   u_INL_2 + u_DNL_2 + u_SGR_2 + u_darkstab_2 + u_noise_2)**0.5
        u_expand = 10 * (self.k * u_1sigma)
        u_ref = np.uint8(np.clip(u_expand, 0, 250))

        return u_ref
    
    def calc_smear_band(self, all_data):
        ### NOT CURRENT IN USE - TO BE UPDATED FOR V1.0###
        # Initialise arrays
        smear_band = np.array([])
        pre_Xsm = np.array([])

        #Loop through bands
        for i, sig in enumerate(all_data):
            #pixel exposure during transfer
            Xsm = (rad_conf.t_trans/rad_conf.n_pxl_tot)*(sig/rad_conf.t_int)

            # Average neighbouring bands per pixel Xsmear values, and use this
            # value to cover the region between the central pixels (at the end
            # of the CCD where this is not possible assume the edge values are
            # 0). Then, add up total signal across the entire CCD.
 
            if i == 0:
                ave_Xsm = (Xsm/2)*((rad_conf.ls[0,0]-390.)/rad_conf.delta_l)
                smear_band = ave_Xsm
                pre_Xsm = Xsm
 
            elif i in range(1, len(rad_conf.ls[0,:])-1):
                ave_Xsm = (Xsm+pre_Xsm)/2*((rad_conf.ls[0,i]-rad_conf.ls[0,i-1])/rad_conf.delta_l)
                smear_band = smear_band + ave_Xsm
                pre_Xsm = Xsm
 
            else:  
                ave_Xsm = (Xsm+pre_Xsm)/2*((1050.-rad_conf.ls[0,-1])/rad_conf.delta_l)
                smear_band = smear_band + ave_Xsm
                pre_Xsm = Xsm

        return smear_band
