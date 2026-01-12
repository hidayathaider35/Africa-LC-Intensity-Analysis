// ======================================================================
// === 1. SETUP: DATA, ROI, AND RECLASSIFICATION (Unchanged) ============
// ======================================================================

// Load your Africa sub-regions shapefile from Assets
// *** RENAME THIS to match your asset path ***
//Please import the IPCC_5_roi_shapefile,  make sure you uploaded it into your assets before importing here
var roi = //IPCC_5_roi_shapefile

// Define the 9 region names from your "Acronym" column
var regionNames = ['EAF', 'MED', 'SAF', 'SAH', 'WAF'];

// Define your reclassification map
// 'From' classes (original values)
var from_classes = [
  10, 11, 12, 20, // 1: Cropland
  51, 52, 61, 62, 71, 72, 81, 82, 91, 92, // 2: Forest
  120, 121, 122, // 3: Shrubland
  130, // 4: Grassland
  140, // 5: Tundra
  181, 182, 183, 184, 185, 186, 187, // 6: Wetland
  190, // 7: Impervious surface
  150, 152, 153, 200, 201, 202, // 8: Bare area
  210, // 9: Waterbody
  220, // 10: Permanent snow and ice
  0    // Mask out 'Filled value'
];

// 'To' classes (your new 1-10 system)
var to_classes = [
  1, 1, 1, 1, // 1: Cropland
  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, // 2: Forest
  3, 3, 3, // 3: Shrubland
  4, // 4: Grassland
  5, // 5: Tundra
  6, 6, 6, 6, 6, 6, 6, // 6: Wetland
  7, // 7: Impervious surface
  8, 8, 8, 8, 8, 8, // 8: Bare area
  9, // 9: Waterbody
  10, // 10: Permanent snow and ice
  0  // Map 0 to 0, to be masked
];

// Reclassification function
function reclassify(image) {
  var remapped = image.remap(from_classes, to_classes).rename('landcover');
  return remapped.updateMask(remapped.neq(0)); // Mask out the 0 values
}

// ======================================================================
// === 2. PREPARE IMAGE COLLECTION (Unchanged) =========================
// ======================================================================

// Load and process the 1985, 1990, 1995 data
var five_year_col = ee.ImageCollection('projects/sat-io/open-datasets/GLC-FCS30D/five-years-map');
var five_year_img = five_year_col.mosaic();
var img1985 = reclassify(five_year_img.select('b1'))
  .set('year', 1985).set('system:time_start', ee.Date('1985-01-01'));
var img1990 = reclassify(five_year_img.select('b2'))
  .set('year', 1990).set('system:time_start', ee.Date('1990-01-01'));
var img1995 = reclassify(five_year_img.select('b3'))
  .set('year', 1995).set('system:time_start', ee.Date('1995-01-01'));
var pre2000_col = ee.ImageCollection.fromImages([img1985, img1990, img1995]);

// Load the annual collection and mosaic it into ONE image with 23 bands (b1-b23)
var annual_mosaic = ee.ImageCollection('projects/sat-io/open-datasets/GLC-FCS30D/annual').mosaic();

// Create a list of years from 2000 to 2022
var years_annual = ee.List.sequence(2000, 2022);

// Map over the list of years to create a new ImageCollection
var annual_col = ee.ImageCollection.fromImages(
  years_annual.map(function(year) {
    year = ee.Number(year);
    var band_index = year.subtract(1999);
    var band_name = ee.String('b').cat(band_index.format('%d'));
    var yearly_image = annual_mosaic.select([band_name]).rename('landcover'); 
    
    return reclassify(yearly_image)
      .set('year', year)
      .set('system:time_start', ee.Date.fromYMD(year, 1, 1));
  })
);

// Merge the two collections
var full_col = pre2000_col.merge(annual_col).sort('system:time_start');

// Define the specific years we will use for our intervals
var yearsToAnalyze = [1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022];

// Filter the full collection to only these years
var analysis_col = full_col.filter(ee.Filter.inList('year', yearsToAnalyze));

Map.addLayer(analysis_col.first(), {min: 1, max: 10, palette: ['#ffff64', '#006400', '#966400', '#ffb432', '#dcdcdc', '#00a884', '#c31400', '#fff5d7', '#0046c8', '#ffffff']}, 'Reclassified 1985');

// ======================================================================
// === 2.5 (NEW) VISUALIZE REGIONS =====================================
// ======================================================================

Map.centerObject(roi, 4); // Center the map on your full ROI

// Create a list of colors for visualization
var colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#00FFFF', '#FF00FF', '#C0C0C0', '#808080', '#800000'];

// Loop through each region name and add it to the map
regionNames.forEach(function(name, index) {
  var region = roi.filter(ee.Filter.eq('LAB', name));
  var color = colors[index % colors.length]; // Get a color
  
  // Add a styled layer to the map
  Map.addLayer(region.style({color: color, fillColor: color + '80'}), // 80 = 50% transparency
    {}, 
    name); // Layer name is the region name
});


// ======================================================================
// === 3. CALCULATE & EXPORT TRANSITION MATRICES (MODIFIED) ===========
// ======================================================================

print('Starting to queue export tasks...');

regionNames.forEach(function(regionName) {
  // Get the geometry for the current region
  var regionGeom = roi.filter(ee.Filter.eq('LAB', regionName)).geometry();
  
  // Iterate through the time intervals
  for (var i = 0; i < yearsToAnalyze.length - 1; i++) {
    var year1 = yearsToAnalyze[i];
    var year2 = yearsToAnalyze[i+1];
    var duration = year2 - year1;
    
    var image1 = analysis_col.filter(ee.Filter.eq('year', year1)).first().clip(regionGeom);
    var image2 = analysis_col.filter(ee.Filter.eq('year', year2)).first().clip(regionGeom);
    
    var transitionImage = image1.multiply(100).add(image2).rename('transition_code');
    var areaImage = ee.Image.pixelArea().addBands(transitionImage);
    
    var reduction = areaImage.reduceRegion({
      reducer: ee.Reducer.sum().group({
        groupField: 1, // 'transition_code' band
        groupName: 'code',
      }),
      geometry: regionGeom,
      scale: 30, // Native resolution of GLC-FCS30D
      maxPixels: 1e13
    });
    
    var transition_list = ee.List(reduction.get('groups'));
    
    var features = transition_list.map(function(item) {
      item = ee.Dictionary(item);
      var code = ee.Number(item.get('code'));
      
      // --- MODIFICATION: Convert m² to km² ---
      var area_m2 = ee.Number(item.get('sum'));
      var area_km2 = area_m2.divide(1000000); // 1,000,000 m² in 1 km²
      // --- END MODIFICATION ---
      
      var from_class = code.divide(100).floor().toInt();
      var to_class = code.mod(100).toInt();
      
      return ee.Feature(null, {
        'region': regionName,
        'year_initial': year1,
        'year_final': year2,
        'interval_duration_yrs': duration,
        'from_class': from_class,
        'to_class': to_class,
        'area_km2': area_km2 // <-- Use the new variable
      });
    });
    
    var transitions_fc = ee.FeatureCollection(features);
    
    var fileName = 'TransitionMatrix_' + regionName + '_' + year1 + '_to_' + year2;
    Export.table.toDrive({
      collection: transitions_fc,
      description: fileName,
      fileNamePrefix: fileName,
      folder:"africa_5r_my_matrix",
      fileFormat: 'CSV'
    });
    
  } // end year loop
}); // end region loop
