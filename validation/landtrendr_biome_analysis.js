// ============================================================
// LandTrendr Forest Disturbance Validation for Avian Pollination MS
// Task-only version — runs as batch job on GEE servers, immune to sleep
// ============================================================

var biomes = ee.FeatureCollection('projects/ee-dzl19960511/assets/biomes');
var lcms = ee.ImageCollection('USFS/GTAC/LCMS/v2021-7');

// Build annual fast-loss area images (m² per pixel)
var annualLoss = ee.ImageCollection(
  ee.List.sequence(1985, 2017).map(function(y) {
    y = ee.Number(y);
    var img = lcms.filterDate(
      ee.Date.fromYMD(y, 1, 1), ee.Date.fromYMD(y, 12, 31)
    ).first();
    return img.select('Change').eq(2)
      .unmask(0)                     // treat CONUS-outside pixels as 0 loss
      .multiply(ee.Image.pixelArea())
      .set('year', y);
  })
);

// reduceRegions: one image → all biomes (avoids concurrency explosion)
var results = annualLoss.map(function(img) {
  var yr = img.get('year');
  return img.reduceRegions({
    collection: biomes,
    reducer: ee.Reducer.sum(),
    scale: 30,
    crs: 'EPSG:5070'
  }).map(function(f) {
    return f.set('year', yr);
  });
}).flatten();

// Export as batch task → runs on GEE servers, survives sleep/network drop
Export.table.toDrive({
  collection: results,
  description: 'LandTrendr_biome_forest_loss_1985_2017',
  fileFormat: 'CSV',
  selectors: ['biome', 'year', 'sum']
});

// Chart preview (uses image-based API, more efficient than FeatureCollection chart)
var chart = ui.Chart.image.seriesByRegion({
  imageCollection: annualLoss,
  regions: biomes,
  reducer: ee.Reducer.sum(),
  scale: 500,              // coarser scale for fast preview
  xProperty: 'year',
  seriesProperty: 'biome'
});
chart.setChartType('LineChart');
chart.setOptions({
  title: 'Forest Fast-Loss Area by Breeding Biome (1985-2017)',
  vAxis: {title: 'Loss area (m²)'},
  hAxis: {title: 'Year', format: '####'},
  lineWidth: 2,
  interpolateNulls: true
});
print(chart);

// Export chart as SVG/PDF to Google Drive
Export.chart.toDrive({
  chart: chart,
  description: 'LandTrendr_biome_forest_loss_chart',
  fileNamePrefix: 'LandTrendr_biome_forest_loss_chart',
  fileFormat: 'SVG',
  dimensions: '1200x600'
});

print('2 tasks submitted. Open the "Tasks" tab (right panel) and click "Run" for both.');
print('The tasks run on Google servers — you can close your laptop.');
