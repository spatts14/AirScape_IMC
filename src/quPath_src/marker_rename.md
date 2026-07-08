# Groovy script for qupath

// QuPath script to rename channels and set pixel calibration to 1 um/px

// ---- 1. Rename channels ----
  def channelNames = [
    "Collagen_I",
    "Pi16",
    "Vimentin",
    "VWF",
    "Collagen_VI",
    "CD16",
    "GATA3",
    "PDGFRb",
    "CD11b",
    "PGP_9.5",
    "CD31",
    "CD45",
    "Collagen_IV",
    "CD11c",
    "CD206",
    "CD4",
    "ECadherin",
    "CD68",
    "SPARC",
    "CD20",
    "CD8a",
    "VEGF",
    "KRT5",
    "PDGFRa",
    "ECP",
    "SCGB1a1",
    "Ki67",
    "MPO",
    "CD3",
    "CD66a",
    "MUC5AC",
    "Tryptase",
    "HLA-DR",
    "Fibronectin",
    "Collagen_III",
    "DNA1",
    "DNA2",
    "aSMA"
  ]
  
  def imageData = getCurrentImageData()
  def server = imageData.getServer()
  def nChannels = server.nChannels()
  
  if (channelNames.size() != nChannels) {
    println "WARNING: Number of channel names (${channelNames.size()}) does not match number of channels in image (${nChannels})"
  }
  
  setChannelNames(*channelNames)
  
  // ---- 2. Set pixel calibration to 1 um x 1 um ----
    setPixelSizeMicrons(1.0, 1.0)
  
  println "Channel names updated successfully!"
  for (int i = 0; i < server.nChannels(); i++) {
    println "Channel ${i+1}: ${server.getChannel(i).getName()}"
  }
  
  println "Pixel calibration set to 1.0 x 1.0 um/px"