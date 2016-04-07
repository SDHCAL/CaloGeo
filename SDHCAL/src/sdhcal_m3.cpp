//====================================================================
//  AIDA Detector description implementation
//  for CALICE SDHCAL m3 prototype
//--------------------------------------------------------------------
//
//  Author     : R. Ete
//
// Basic idea:
// 1. Create the Hcal module envelope.
//    Note: with default material air.
//    
// 2. Create the Hcal Layer.
//    Create the Layer with slices (Steel235,Polystyrene,Cu,FR4,air).
//    Place each slice into the layer with the right position,
//    And registry the ID for slices
//
// 3. Place the same Layer into the module envelope.
//    It will be repeated, and registry the IDs for layer.
//
// 4. Place the module into the world volume,
//    with the right position and rotation.
//    And registry the ID for System.
//
// 5. Customer material FR4 and Steel235 defined in materials.xml
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"

#include "DDSegmentation/CartesianGridXY.h"
#include "DDSegmentation/Segmentation.h"
#include "DDSegmentation/SegmentationFactory.h"

using namespace DD4hep;
using namespace DD4hep::Geometry;
using namespace DD4hep::DDSegmentation;

class SDHCALSegmentation : public DDSegmentation::CartesianGridXY
{
public:
	/// Default constructor passing the encoding string
	SDHCALSegmentation(const std::string &cellIEncoding = "");
	/// destructor
	~SDHCALSegmentation();

	/// convert cellID to position
	DDSegmentation::Vector3D position(const DDSegmentation::CellID& cID) const;

	///
	DDSegmentation::CellID cellID(const DDSegmentation::Vector3D& localPosition, const DDSegmentation::Vector3D& globalPosition, const DDSegmentation::VolumeID& vID) const;

private:
	std::string _difId;
	std::string _asicId;
	std::string _channelId;
};

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

inline SDHCALSegmentation::SDHCALSegmentation(const std::string &cellEncoding) :
		DDSegmentation::CartesianGridXY(cellEncoding)
{
	_type = "SDHCALSegmentation";
	_description = "SDHCAL specific segmentation given a layer ID, a dif ID, an asic ID and channel ID";

	registerIdentifier("identifier_dif", "Cell ID identifier for dif id", _difId, "Dif_Id");
	registerIdentifier("identifier_asic", "Cell ID identifier for asic id", _asicId, "Asic_Id");
	registerIdentifier("identifier_channel", "Cell ID identifier for channel id", _channelId, "Channel_Id");
}

//--------------------------------------------------------------------------------------------

inline SDHCALSegmentation::~SDHCALSegmentation()
{
}

//--------------------------------------------------------------------------------------------

inline DDSegmentation::Vector3D SDHCALSegmentation::position(const DDSegmentation::CellID& cID) const
{
	// TODO decode position as a function of the cellID for sdhcal
	return DDSegmentation::CartesianGridXY::position(cID);
}

//--------------------------------------------------------------------------------------------

inline DDSegmentation::CellID SDHCALSegmentation::cellID(const DDSegmentation::Vector3D& localPosition,
		const DDSegmentation::Vector3D& globalPosition, const DDSegmentation::VolumeID& vID) const
{
	// TODO decode cellID as a function of the positions for sdhcal
	return DDSegmentation::CartesianGridXY::cellID(localPosition, globalPosition, vID);
}

//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

static Ref_t create_detector(LCDD& lcdd, xml_h element, SensitiveDetector sens)
{
  xml_det_t   x_det     = element;
  Layering    layering(x_det.child(_U(layers)));
  std::string det_name    = x_det.nameStr();
  std::string det_type    = x_det.typeStr();
  Material    air         = lcdd.air();
  xml_dim_t   dim         = x_det.dimensions();
  xml_dim_t   x_pos  (x_det.child(_U(position)));
  xml_dim_t   x_rot  (x_det.child(_U(rotation)));

  Translation3D det_pos(x_pos.x(),x_pos.y(),x_pos.z());
  //Rotation is  ROOT::Math::RotationZYX
  //The input is Rotation(z,y,x)
  RotationZYX   det_rot(x_rot.z(),x_rot.y(),x_rot.x());
  
  DetElement   sdet(det_name,x_det.id());
  Volume      motherVol = lcdd.pickMotherVolume(sdet); 


  
  // ========= Create Hcal Modules envelope ============================
  //  They will be the volume for placing the Hcal Layers.
  //  Themselves will be placed into the world volume.
  // ==========================================================================
  
  // Hcal module shape
  double box_half_x  = dim.x()/2.0;
  double box_half_y  = dim.y()/2.0;
  double box_half_z  = dim.z()/2.0;
 
  Box    BoxModule(box_half_x,box_half_y,box_half_z);
  
  // define the name of  Module
  std::string envelopeVol_name   = det_name+_toString("_envelope");
  
  Volume envelopeVol(envelopeVol_name,BoxModule,air);
  
  // Set envelope volume attributes.
  envelopeVol.setAttributes(lcdd,x_det.regionStr(),x_det.limitsStr(),x_det.visStr());
  
  
  // ========= Create Hcal Layers===== ==============================
  // It will be the sub volume for placing the slices.
  // Itself will be placed into the Hcal modules envelope.
  // ================================================================
  
  // create Layer (air) and place the slices into it. 
  // place the Layer into the Hcal Modules envelope.
  
  // Hcal layer start position
  double layer_pos_z     = - box_half_z;                      
  
  // Create Hcal Chamber without radiator
  // Place into the Hcal module envelope 
  int layer_num = 0;
  int module_num = 0;

  for(xml_coll_t c(x_det.child(_U(layers)),_U(layer)); c; ++c)
  {
    xml_comp_t   x_layer = c;

    //Layering    layering(x_layer);
    int          repeat = x_layer.repeat();
    const Layer* lay    = layering.layer(layer_num); // Get the layer from the layering engine.
    
    std::string layer_name      = det_name+ _toString(module_num,"_module%d_layer");
    double layer_thickness = lay->thickness();
    DetElement  layer(layer_name,"layerModule",x_det.id());
    
    // Layer box & volume
    double layer_dim_x = box_half_x;
    double layer_dim_y = box_half_y;
    double layer_dim_z = layer_thickness/2.0;
    
    // Build chamber including air gap
    // The Layer will be filled with slices, 
    Volume layer_vol(layer_name, Box(layer_dim_x,layer_dim_y,layer_dim_z), air);
    
    
    // ========= Create sublayer slices =========================================
    // Create and place the slices into Layer
    // ==========================================================================
    
    // Create the slices (sublayers) within the Hcal Chamber.
    double slice_pos_z = -(layer_thickness / 2.0);
    int slice_number = 0;

    for(xml_coll_t k(x_layer,_U(slice)); k; ++k)
    {
      xml_comp_t x_slice = k;
      std::string   slice_name      = layer_name + _toString(slice_number,"_slice%d");
      double   slice_thickness = x_slice.thickness();
      Material slice_material  = lcdd.material(x_slice.materialStr());

      DetElement slice(layer,_toString(slice_number,"slice%d"),x_det.id());
      
      slice_pos_z += slice_thickness / 2.0;
      
      // Slice volume & box
      Volume slice_vol(slice_name,Box(layer_dim_x,layer_dim_y,slice_thickness/2.0),slice_material);
      
      if ( x_slice.isSensitive() )
      {
		sens.setType("calorimeter");
		slice_vol.setSensitiveDetector(sens);
      }

      // Set region, limitset, and vis.
      slice_vol.setAttributes(lcdd,x_slice.regionStr(),x_slice.limitsStr(),x_slice.visStr());
      // slice PlacedVolume
      PlacedVolume slice_phv = layer_vol.placeVolume(slice_vol,Position(0,0,slice_pos_z));
      slice_phv.addPhysVolID("slice",slice_number);
      
      slice.setPlacement(slice_phv);
      // Increment Z position for next slice.
      slice_pos_z += slice_thickness / 2.0;
      // Increment slice number.
      ++slice_number;             
    }
    // Set region, limitset, and vis.
    layer_vol.setAttributes(lcdd,x_layer.regionStr(),x_layer.limitsStr(),x_layer.visStr());
    
    
    // ========= Place the Layer (i.e. Chamber) =================================
    // Place the Layer into the Hcal module envelope.
    // with the right position, and registry the ID layer
    // ==========================================================================
    
    for (int j = 0; j < repeat; j++)
    {
      
      // Layer position in z within the Modules.
      layer_pos_z += layer_thickness / 2.0;
      
      PlacedVolume layer_phv = envelopeVol.placeVolume(layer_vol,Position(0,0,layer_pos_z));

      // registry the ID of Layer
      layer_phv.addPhysVolID("layer",layer_num);
      layer_phv.addPhysVolID("module",j);

      // then setPlacement for it.
      layer.setPlacement(layer_phv);
      
      // Increment the layer_pos_z
      layer_pos_z += layer_thickness / 2.0;
      ++layer_num;         
    }
    
    module_num++; 
  }
  
  
  // for the alignment in the compact XML file
  Transform3D Tr3D = Transform3D(det_rot,det_pos);
  PlacedVolume env_phv = motherVol.placeVolume(envelopeVol,Tr3D);
  
  
  env_phv.addPhysVolID("system",x_det.id());
  sdet.setPlacement(env_phv);
  
  return sdet;
}

// register sdhcal segmentation
REGISTER_SEGMENTATION(SDHCALSegmentation)

// register sdhcal detector m3 prototype
DECLARE_DETELEMENT(SDHcal_detector, create_detector)
