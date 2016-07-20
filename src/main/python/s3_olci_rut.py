# -*- coding: utf-8 -*-
"""
Created on Thurs May 26 17:51:33 2016

@author: shunt
"""

import snappy
import s3_olci_rut_algo
import numpy as np
import datetime
import s3_olci_l1_rad_conf as rad_conf

S3_OLCI_TYPE_STRING = 'S3_OLCI_Level-1B'

class S3OLCIRutOp:
    def __init__(self):
        self.source_product = None
        self.product_meta = None
        self.datastrip_meta = None
        self.rut_algo = s3_olci_rut_algo.S3OLCIRutAlgo()
        self.unc_band = None
        self.time_init = datetime.datetime(2016, 2, 16, 10, 00)  # S3A launch date 16-feb-2016, time is indifferent
        self.sourceBandMap = None

    def initialize(self, context):
        self.source_product = context.getSourceProduct()

        if self.source_product.getProductType() != S3_OLCI_TYPE_STRING:
            raise RuntimeError('Source product must be of type "' + S3_OLCI_TYPE_STRING + '"')

        self.product_meta, self.datastrip_meta, granules_meta = self.source_product.getMetadataRoot().getElements()

        # todo - check if there is a granule

        self.toa_band_names = context.getParameter('band_names')
        self.rut_algo.k = self.get_k(context)
        self.rut_algo.unc_select = self.get_unc_select(context)

        scene_width = self.source_product.getSceneRasterWidth()
        scene_height = self.source_product.getSceneRasterHeight()

        rut_product = snappy.Product(self.source_product.getName() + '_rut', 'S3_OLCI_RUT', scene_width, scene_height)
        snappy.ProductUtils.copyGeoCoding(self.source_product, rut_product)
        self.sourceBandMap = {}
        for name in self.toa_band_names:
            source_band = self.source_product.getBand(name)
            unc_toa_band = snappy.Band(name + '_rut', snappy.ProductData.TYPE_UINT8, source_band.getRasterWidth(),
                                       source_band.getRasterHeight())
            unc_toa_band.setDescription('Uncertainty of ' + name + ' (coverage factor k=' + str(self.rut_algo.k) + ')')
            unc_toa_band.setNoDataValue(250)
            unc_toa_band.setNoDataValueUsed(True)
            rut_product.addBand(unc_toa_band)
            self.sourceBandMap[unc_toa_band] = source_band

        context.setTargetProduct(rut_product)

    def computeTile(self, context, band, tile):
        source_band = self.sourceBandMap[band]
        toa_band_id = source_band.getSpectralBandIndex() - 1

        toa_tile = context.getSourceTile(source_band, tile.getRectangle())
        toa_samples = toa_tile.getSamplesFloat()

        # this is the core where the uncertainty calculation should grow
        unc = self.rut_algo.unc_calculation(np.array(toa_samples, dtype=np.uint16), toa_band_id)

        tile.setSamples(unc)

    def get_k(self, context):
        return (context.getParameter('coverage_factor'))

    def get_unc_select(self, context):
        return ([context.getParameter('Instrument_noise'), context.getParameter('Instrument_aging'),
                 context.getParameter('CCD_stability'), context.getParameter('Period_signal_error'),
                 context.getParameter('Diffuser_characterisation'), context.getParameter('Diffuser_modelisation'),
                 context.getParameter('Diffuser_alignment'), context.getParameter('Calibration_diffuser_straylight'),
                 context.getParameter('Calibration_camera_straylight'), context.getParameter('Calibration_diffuser_aging'),
                 context.getParameter('Reference_diffuser_aging'), context.getParameter('Video_chain_non-linearity'),
                 context.getParameter('ADC_non-linearity'), context.getParameter('Offset_compensation'),
                 context.getParameter('Dark_stability'), context.getParameter('Smear_gain_compensation')])
