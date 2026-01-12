// ======================================================================
// === SPATIALIZING INTENSITY ANALYSIS: AFRICA (1985-2022) =============
// ======================================================================

// 1. SETUP
// *** RENAME 'roi' to match your asset variable name ***
var roi = //IPCC_5_roi_shapefile 

// Define reclassification (Original 1-10 Scheme)
// 1:CRP, 2:FST, 3:SHR, 4:GRS, 5:TUD, 6:WET, 7:IMP, 8:BAL, 9:WTR, 10:PSI
var from_classes = [10, 11, 12, 20, 51, 52, 61, 62, 71, 72, 81, 82, 91, 92, 120, 121, 122, 130, 140, 181, 182, 183, 184, 185, 186, 187, 190, 150, 152, 153, 200, 201, 202, 210, 220, 0];
var to_classes =   [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 7, 8, 8, 8, 8, 8, 8, 9, 10, 0];

function reclassify(image) {
  return image.remap(from_classes, to_classes).rename('landcover').updateMask(image.neq(0));
}

// 2. LOAD DATA (1985 Start vs 2022 End)
var five_year_col = ee.ImageCollection('projects/sat-io/open-datasets/GLC-FCS30D/five-years-map');
var img1985 = reclassify(five_year_col.mosaic().select('b1')).clip(roi);

var annual_mosaic = ee.ImageCollection('projects/sat-io/open-datasets/GLC-FCS30D/annual').mosaic();
var img2022 = reclassify(annual_mosaic.select('b23')).clip(roi);

// 3. CREATE TRANSITION MAP
// Formula: Initial * 100 + Final. 
var transition = img1985.multiply(100).add(img2022).rename('transition');

// ======================================================================
// 4. MAP THE EXPANDED TRANSITIONS
// ======================================================================

// Define the remapping for visualization
var target_transitions = transition.remap(
  // From-To Codes (Initial*100 + Final):
  [
    201,  // FST to CRP: Deforestation
    
    203,  // FST to SHR: Forest Degradation 
    
    301, 401, 801, // SHR/GRS/BAL to CRP: Ag Expansion
    
    206,  // FST to WET: Wetland Expansion
    
    304,  // SHR to GRS: Shrub Degradation
    
    408   // GRS to BAL: Desertification (NEW ADDITION)
  ],
  // Visualization IDs (Legend Categories):
  [
    1,    // Deforestation
    
    2,    // Forest Degradation
    
    3, 3, 3, // Ag Expansion
    
    4,    // Wetland Expansion
    
    5,    // Shrub Degradation
    
    6     // Desertification
  ],
  0 // Default value (masked)
);

var map_layer = target_transitions.updateMask(target_transitions.neq(0));

// ======================================================================
// 5. CALCULATE STATISTICS FOR CSV EXPORT
// ======================================================================

// Calculate area in km2 for each class (1-6)
var areaImage = ee.Image.pixelArea().divide(1e6).addBands(map_layer);

var stats = areaImage.reduceRegion({
  reducer: ee.Reducer.sum().group({
    groupField: 1,
    groupName: 'class_id',
  }),
  geometry: roi.geometry(),
  scale: 1000, // Resolution for calculation (adjust to 30 for high precision, 1000 for speed)
  maxPixels: 1e13
});

// Create Features for CSV Export
var statList = ee.List(stats.get('groups'));
var classNames = ee.Dictionary({
  '1': 'Deforestation (FST->CRP)',
  '2': 'Forest Degradation (FST->SHR)',
  '3': 'Ag. Expansion (SHR/GRS/BAL->CRP)',
  '4': 'Wetland Expansion (FST->WET)',
  '5': 'Shrub Degradation (SHR->GRS)',
  '6': 'Desertification (GRS->BAL)'
});

var exportFeatures = statList.map(function(item) {
  var dict = ee.Dictionary(item);
  // FIX: Cast to ee.Number before formatting
  var classId = ee.Number(dict.get('class_id')).format('%d'); 
  var area = dict.get('sum');
  var name = classNames.get(classId);
  return ee.Feature(null, {
    'Class_ID': classId,
    'Transition_Type': name,
    'Area_km2': area,
    'Region': 'AFRICA_Overall',
    'Period': '1985-2022'
  });
});

var exportFC = ee.FeatureCollection(exportFeatures);

// Export CSV Task
Export.table.toDrive({
  collection: exportFC,
  description: 'Africa_Transition_Statistics_1985_2022',
  folder: 'Africa_Intensity_Analysis',
  fileFormat: 'CSV',
  selectors: ['Transition_Type', 'Area_km2', 'Class_ID', 'Region', 'Period']
});

// ======================================================================
// 6. VISUALIZATION
// ======================================================================
Map.centerObject(roi, 4);

// Background: Grey-scale version of 2022 land cover for context
var greyVis = {
  min: 1, 
  max: 10, 
  palette: ['e0e0e0', 'd0d0d0', 'c0c0c0', 'b0b0b0', 'a0a0a0', '909090', '808080', '707070', '606060', 'ffffff']
};
Map.addLayer(img2022, greyVis, 'Background Context (2022)', true, 0.6);

// The Transition Colors
var transVis = {
  min: 1, 
  max: 6, 
  palette: [
    '#FF0000', // 1: Deforestation (Red)
    '#FFA500', // 2: Forest Degradation (Orange)
    '#FFFF00', // 3: Ag Expansion (Yellow)
    '#00CED1', // 4: Wetland Expansion (Turquoise)
    '#8B4513', // 5: Shrub Degradation (Brown)
    '#F0E68C'  // 6: Desertification (Khaki/Sand)
  ]
};

Map.addLayer(map_layer, transVis, 'Major Transitions (1985-2022)');

// ======================================================================
// 7. LEGEND (With Areas)
// ======================================================================
// Note: Legend values are derived from previous analysis summation.
// Run the Export Task to get the exact pixel-count precision.

var legend = ui.Panel({
  style: {
    position: 'bottom-left', 
    padding: '8px 15px',
    backgroundColor: 'white'
  }
});

var title = ui.Label({
  value: 'Major Transitions (1985-2022)', 
  style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 8px 0'}
});
legend.add(title);

var subtitle = ui.Label({
  value: 'Areas are approx. totals', 
  style: {fontSize: '10px', color: 'gray', margin: '0 0 8px 0'}
});
legend.add(subtitle);

var addLegendRow = function(color, label, area) {
  var box = ui.Label({
    style: {
      backgroundColor: color, 
      padding: '8px', 
      margin: '0 4px 4px 0',
      border: '1px solid black'
    }
  });
  var desc = ui.Label({
    value: label + ' (' + area + ')', 
    style: {margin: '0', fontSize: '11px'}
  });
  return ui.Panel({
    widgets: [box, desc], 
    layout: ui.Panel.Layout.Flow('horizontal')
  });
};

legend.add(addLegendRow('#FF0000', 'Deforestation (FST→CRP)', '400K km²'));
legend.add(addLegendRow('#FFA500', 'Forest Degradation (FST→SHR)', '1.13M km²'));
legend.add(addLegendRow('#FFFF00', 'Ag. Expansion (All→CRP)', '792K km²'));
legend.add(addLegendRow('#00CED1', 'Wetland Expansion (FST→WET)', '285K km²'));
legend.add(addLegendRow('#8B4513', 'Shrub Degradation (SHR→GRS)', '683K km²'));
legend.add(addLegendRow('#F0E68C', 'Desertification (GRS→BAL)', '317K km²'));

Map.add(legend);

// 8. EXPORT MAP IMAGE
Export.image.toDrive({
  image: map_layer,
  description: 'Africa_Transitions_Map_1985_2022_v2',
  folder: 'Africa_Intensity_Analysis',
  region: roi.geometry(),
  scale: 1000, 
  crs: 'EPSG:4326',
  maxPixels: 1e13
});
