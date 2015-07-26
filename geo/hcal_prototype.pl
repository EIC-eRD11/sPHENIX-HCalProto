#############################################################
#############################################################
#  Build the sPHENIX HCAL prototype detector
#  All dimensions are in mm
#############################################################

use strict;
use warnings;
use Getopt::Long;
use Math::Trig;

our %configuration;

my $DetectorMother="root";
my $DetectorName = 'rhic_sphnx_hcalproto';

# no overlap gap
my $no_overlap = 0.00001;

#######------ Define the holder Box for Detectors ------#######
my $box_halfx = 500;
my $box_halfy = 700;
my $box_halfz = 1400;

my $offset = $box_halfz;

my $box_name = "box_holder";
my $box_mat = "G4_AIR";
sub build_box()
{
	my @box_pos  = ( 0.0, 0.0, $offset );
	my @box_size = ( $box_halfx, $box_halfy, $box_halfz );
        my %detector=init_det();
        $detector{"name"} = "$DetectorName\_$box_name";
        $detector{"mother"} = "$DetectorMother";
        $detector{"description"} = "$DetectorName\_$box_name";
        $detector{"pos"} = "$box_pos[0]*mm $box_pos[1]*mm $box_pos[2]*mm";
        $detector{"color"} = "ffffff";
        $detector{"type"} = "Box";
        $detector{"visible"} = "1";
	$detector{"rotation"} = "0*deg 0*deg 0*deg";
	$detector{"dimensions"} = "$box_size[0]*mm $box_size[1]*mm $box_size[2]*mm";
        $detector{"material"} = "$box_mat";
        $detector{"sensitivity"} = "no";
        $detector{"hit_type"}    = "no";
        $detector{"identifiers"} = "no";
        print_det(\%configuration, \%detector);
}


## delta theta angle for each cell (scint tile + steel plate) 
my $in_delta_theta = 1.25*pi/180.; ## tune this para to be a smaller one
# minmum theta angle
my $in_theta_min = -12.5*pi/180.;

# 1164.5 is the inner radius, 1370 is outer radius
my @hcal_in_radi = (1164.5, 1370.0);
my $hcal_in_radimid = ($hcal_in_radi[0]+$hcal_in_radi[1])/2.;
my $hcal_in_deep = $hcal_in_radi[1]-$hcal_in_radi[0];
my $in_tile_angle = 32.*pi/180.;
# 0.5*$hcal_in_deep/cos($in_tile_angle) in Z
# 0.5*14.5 (steel) + 0.5*8 (scint) in Y
# 247.5+1 in X (outer Z), 
my @hcal_in_halfbox = (248.5, 11.25, 0.5*202.89+1); 

my $in_abs_pDz = 0.5*202.89;
my $in_abs_pTheta = 0.563;  ### calculated with absorber dims
my $in_abs_pPhi = 0;
my $in_abs_pDy1 = 0.5*12.59*sin(54.8*pi/180.);
my $in_abs_pDx1 = $hcal_in_halfbox[0]-1;
my $in_abs_pDx2 = $hcal_in_halfbox[0]-1;
my $in_abs_pAlp1 = 0;
my $in_abs_pAlp2 = 0;
my $in_abs_pDy2 = 0.5*14.27;
my $in_abs_pDx3 = $hcal_in_halfbox[0]-1;
my $in_abs_pDx4 = $hcal_in_halfbox[0]-1;


my $in_tile1_pDz = 0.5*198.1;
my $in_tile1_pTheta = 2.23;  ### calculated with tile1 dims
my $in_tile1_pPhi = 0;
my $in_tile1_pDy1 = 0.5*7.0;
my $in_tile1_pDx1 = 0.5*105.9;
my $in_tile1_pDx2 = 0.5*105.9;
my $in_tile1_pAlp1 = 0;
my $in_tile1_pAlp2 = 0;
my $in_tile1_pDy2 = 0.5*7.0;
my $in_tile1_pDx3 = 0.5*121.3;
my $in_tile1_pDx4 = 0.5*121.3;

my $in_tile2_pDz = 0.5*198.1; 
my $in_tile2_pTheta = 6.67; ### calculated with tile1 and tile2 dims
my $in_tile2_pPhi = 0;
my $in_tile2_pDy1 = 0.5*7.0;
my $in_tile2_pDx1 = 0.5*110.59;
my $in_tile2_pDx2 = 0.5*110.59;
my $in_tile2_pAlp1 = 0;
my $in_tile2_pAlp2 = 0;
my $in_tile2_pDy2 = 0.5*7.0;
my $in_tile2_pDx3 = 0.5*(141.5-15.4);
my $in_tile2_pDx4 = 0.5*(141.5-15.4);

#my $in_abs_mat = "G4_SS310";   # not defined
my $in_abs_mat = "G4_Cu";
my $in_scint_mat = "ScintillatorB";
sub build_hcal_in()
{
	#### contruct holder boxes for each cell
	my %detector;
	for (my $inbox=0; $inbox<20; $inbox++){
		my $inbox_posX = 0.;
		my $inbox_posY = $hcal_in_radimid*sin($in_theta_min+$inbox*$in_delta_theta+$in_delta_theta/2.);
		my $inbox_posZ = $hcal_in_radimid*cos($in_theta_min+$inbox*$in_delta_theta+$in_delta_theta/2.) - $offset;
		%detector=init_det();
                $detector{"name"} = "$DetectorName\_$box_name\_HCalIn_$inbox";
                $detector{"mother"} = "$DetectorName\_$box_name";
                $detector{"description"} = "$DetectorName\_$box_name\_HCalIn_$inbox";
                $detector{"pos"} = "$inbox_posX*mm $inbox_posY*mm $inbox_posZ*mm";
                $detector{"color"} = "ffffff";
                $detector{"type"} = "Box";
                $detector{"dimensions"} = "$hcal_in_halfbox[0]*mm $hcal_in_halfbox[1]*mm $hcal_in_halfbox[2]*mm";
                $detector{"material"} = $box_mat;
                $detector{"rotation"} = "-32*deg 0*deg 0*deg";
                $detector{"visible"} = "1";
                $detector{"sensitivity"} = "no";
                $detector{"hit_type"}    = "no";
                $detector{"identifiers"} = "no";
                print_det(\%configuration, \%detector);

		#### contruct one Cu absorber plate in the holder box
		my $inbox_abs_posX = 0.0;
		my $inbox_abs_posY = -0.5*22.27 + 0.5*($in_abs_pDy1 + $in_abs_pDy2) + $no_overlap;
		my $inbox_abs_posZ = 0.0;
		%detector=init_det();
                $detector{"name"} = "$DetectorName\_$box_name\_HCalIn_$inbox\_Abs";
                $detector{"mother"} = "$DetectorName\_$box_name\_HCalIn_$inbox";
                $detector{"description"} = "$DetectorName\_$box_name\_HCalIn_$inbox\_Abs";
                $detector{"pos"} = "$inbox_abs_posX*mm $inbox_abs_posY*mm $inbox_abs_posZ*mm";
                $detector{"color"} = "ffff00";
                $detector{"type"} = "G4Trap";
                $detector{"dimensions"} = "$in_abs_pDz*mm $in_abs_pTheta*deg $in_abs_pPhi*deg $in_abs_pDy1*mm $in_abs_pDx1*mm $in_abs_pDx2*mm $in_abs_pAlp1*deg $in_abs_pDy2*mm $in_abs_pDx3*mm $in_abs_pDx4*mm $in_abs_pAlp2*mm";
                $detector{"material"} = $in_abs_mat;
		$detector{"rotation"} = "0*deg 0*deg 0*deg";
		#$detector{"visible"} = "1";
		$detector{"sensitivity"} = "no";
		$detector{"hit_type"}    = "no";
		$detector{"identifiers"} = "no";
		print_det(\%configuration, \%detector);

		#### contruct four scint sheets in the holder box
		my @inbox_scint_posX = (56.8, -56.8, 169.0+3.9, -169.0-3.9);
		my $inbox_scint_posY = 0.5*22.27 - 0.5*8; ## half scint width
		my $inbox_scint_posZ = 0.0;
		for (my $insheet=0; $insheet<4; $insheet++){
			%detector=init_det();
			$detector{"name"} = "$DetectorName\_$box_name\_HCalIn_$inbox\_Scnt\_$insheet";
			$detector{"mother"} = "$DetectorName\_$box_name\_HCalIn_$inbox";
			$detector{"description"} = "$DetectorName\_$box_name\_HCalIn_$inbox\_Scnt\_$insheet";
			$detector{"pos"} = "$inbox_scint_posX[$insheet]*mm $inbox_scint_posY*mm $inbox_scint_posZ*mm";
			$detector{"color"} = "00ff00";
			$detector{"type"} = "G4Trap";
			if($insheet==0 || $insheet==1){
				$detector{"dimensions"} = "$in_tile1_pDz*mm $in_tile1_pTheta*deg $in_tile1_pPhi*deg $in_tile1_pDy1*mm $in_tile1_pDx1*mm $in_tile1_pDx2*mm $in_tile1_pAlp1*deg $in_tile1_pDy2*mm $in_tile1_pDx3*mm $in_tile1_pDx4*mm $in_tile1_pAlp2*mm";
			}
			elsif($insheet==2 || $insheet==3){
				$detector{"dimensions"} = "$in_tile2_pDz*mm $in_tile2_pTheta*deg -$in_tile2_pPhi*deg $in_tile2_pDy1*mm $in_tile2_pDx1*mm $in_tile2_pDx2*mm $in_tile2_pAlp1*deg $in_tile2_pDy2*mm $in_tile2_pDx3*mm $in_tile2_pDx4*mm $in_tile2_pAlp2*mm";
			}
			$detector{"material"} = $in_scint_mat;
			if($insheet==1 || $insheet==3){ $detector{"rotation"} = "0*deg 0*deg 180*deg";}
			elsif($insheet==0 || $insheet==2){ $detector{"rotation"} = "0*deg 0*deg 0*deg";}
			$detector{"visible"} = "1";
			$detector{"sensitivity"} = "no";
			$detector{"hit_type"}    = "flux";
			$detector{"identifiers"} = "no";
			print_det(\%configuration, \%detector);
		}

	}
}


## delta theta angle for each cell (scint tile + steel plate) 
my $out_delta_theta = 1.35*pi/180.; ## tune this para to be a smaller one
# minmum theta angle
my $out_theta_min = -13.5*pi/180.;

# 1830 is the inner radius, 2685 is outer radius
my @hcal_out_radi = (1830, 2685);
my $hcal_out_radimid = ($hcal_out_radi[0]+$hcal_out_radi[1])/2.;
my $hcal_out_deep = $hcal_out_radi[1]-$hcal_out_radi[0];
my $out_tile_angle = 12.*pi/180.;
my @hcal_out_halfbox = (488.5, 25.5, 0.5*831.+1); 

my $out_abs_pDz = 0.5*(27.1*cos(74.72*pi/180.)+823.);
my $out_abs_pTheta = 0.563;  ### calculated with absorber dims
my $out_abs_pPhi = 0;
my $out_abs_pDy1 = 0.5*27.1*sin(74.72*pi/180.);
my $out_abs_pDx1 = $hcal_out_halfbox[0]-1;
my $out_abs_pDx2 = $hcal_out_halfbox[0]-1;
my $out_abs_pAlp1 = 0;
my $out_abs_pAlp2 = 0;
my $out_abs_pDy2 = 0.5*42.5;
my $out_abs_pDx3 = $hcal_out_halfbox[0]-1;
my $out_abs_pDx4 = $hcal_out_halfbox[0]-1;

my $out_tile1_pDz = 0.5*828.9;
my $out_tile1_pTheta = 2.57;  ### calculated with tile1 dims
my $out_tile1_pPhi = 0;
my $out_tile1_pDy1 = 0.5*7.0;
my $out_tile1_pDx1 = 0.5*166.25;
my $out_tile1_pDx2 = 0.5*166.25;
my $out_tile1_pAlp1 = 0;
my $out_tile1_pAlp2 = 0;
my $out_tile1_pDy2 = 0.5*7.0;
my $out_tile1_pDx3 = 0.5*240.54;
my $out_tile1_pDx4 = 0.5*240.54;

my $out_tile2_pDz = 0.5*828.9; 
my $out_tile2_pTheta = 7.69; ### calculated with tile1 and tile2 dims
my $out_tile2_pPhi = 0;
my $out_tile2_pDy1 = 0.5*7.0;
my $out_tile2_pDx1 = 0.5*171.5;
my $out_tile2_pDx2 = 0.5*171.5;
my $out_tile2_pAlp1 = 0;
my $out_tile2_pAlp2 = 0;
my $out_tile2_pDy2 = 0.5*7.0;
my $out_tile2_pDx3 = 0.5*(320.95-74.3);
my $out_tile2_pDx4 = 0.5*(320.95-74.3);

my $out_abs_mat = "G4_Fe";
my $out_scint_mat = "ScintillatorB";

sub build_hcal_out()
{
	my %detector;
	for (my $outbox=0; $outbox<20; $outbox++){
		my $outbox_posX = 0.;
		my $outbox_posY = $hcal_out_radimid*sin($out_theta_min+$outbox*$out_delta_theta+$out_delta_theta/2.);
		my $outbox_posZ = $hcal_out_radimid*cos($out_theta_min+$outbox*$out_delta_theta+$out_delta_theta/2.) - $offset;
		%detector=init_det();
                $detector{"name"} = "$DetectorName\_$box_name\_HCalOut_$outbox";
                $detector{"mother"} = "$DetectorName\_$box_name";
                $detector{"description"} = "$DetectorName\_$box_name\_HCalOut_$outbox";
                $detector{"pos"} = "$outbox_posX*mm $outbox_posY*mm $outbox_posZ*mm";
                $detector{"color"} = "ffffff";
                $detector{"type"} = "Box";
                $detector{"dimensions"} = "$hcal_out_halfbox[0]*mm $hcal_out_halfbox[1]*mm $hcal_out_halfbox[2]*mm";
                $detector{"material"} = $box_mat;
                $detector{"rotation"} = "12*deg 0*deg 0*deg";
                $detector{"visible"} = "1";
                $detector{"sensitivity"} = "no";
                $detector{"hit_type"}    = "no";
                $detector{"identifiers"} = "no";
                print_det(\%configuration, \%detector);

		#### contruct one Fe absorber plate in the holder box
		my $outbox_abs_posX = 0.0;
		my $outbox_abs_posY = -0.5*(42.5+8) + 0.5*($out_abs_pDy1 + $out_abs_pDy2) + $no_overlap;
		my $outbox_abs_posZ = 0.0;
		%detector=init_det();
                $detector{"name"} = "$DetectorName\_$box_name\_HCalOut_$outbox\_Abs";
                $detector{"mother"} = "$DetectorName\_$box_name\_HCalOut_$outbox";
                $detector{"description"} = "$DetectorName\_$box_name\_HCalOut_$outbox\_Abs";
                $detector{"pos"} = "$outbox_abs_posX*mm $outbox_abs_posY*mm $outbox_abs_posZ*mm";
                $detector{"color"} = "ffff00";
                $detector{"type"} = "G4Trap";
                $detector{"dimensions"} = "$out_abs_pDz*mm $out_abs_pTheta*deg $out_abs_pPhi*deg $out_abs_pDy1*mm $out_abs_pDx1*mm $out_abs_pDx2*mm $out_abs_pAlp1*deg $out_abs_pDy2*mm $out_abs_pDx3*mm $out_abs_pDx4*mm $out_abs_pAlp2*mm";
                $detector{"material"} = $out_abs_mat;
		$detector{"rotation"} = "0*deg 0*deg 0*deg";
		$detector{"visible"} = "1";
		$detector{"sensitivity"} = "no";
		$detector{"hit_type"}    = "no";
		$detector{"identifiers"} = "no";
		print_det(\%configuration, \%detector);

		#### contruct four scint sheets in the holder box
		my @outbox_scint_posX = (101.7, -101.7, 289.4+18.6, -289.4-18.6);
		my $outbox_scint_posY = 0.5*(42.5+8) - 0.5*8; ## half scint width
		my $outbox_scint_posZ = 0.0;
		for (my $outsheet=0; $outsheet<4; $outsheet++){
			%detector=init_det();
			$detector{"name"} = "$DetectorName\_$box_name\_HCalOut_$outbox\_Scnt\_$outsheet";
			$detector{"mother"} = "$DetectorName\_$box_name\_HCalOut_$outbox";
			$detector{"description"} = "$DetectorName\_$box_name\_HCalOut_$outbox\_Scnt\_$outsheet";
			$detector{"pos"} = "$outbox_scint_posX[$outsheet]*mm $outbox_scint_posY*mm $outbox_scint_posZ*mm";
			$detector{"color"} = "00ff00";
			$detector{"type"} = "G4Trap";
			if($outsheet==0 || $outsheet==1){
				$detector{"dimensions"} = "$out_tile1_pDz*mm $out_tile1_pTheta*deg $out_tile1_pPhi*deg $out_tile1_pDy1*mm $out_tile1_pDx1*mm $out_tile1_pDx2*mm $out_tile1_pAlp1*deg $out_tile1_pDy2*mm $out_tile1_pDx3*mm $out_tile1_pDx4*mm $out_tile1_pAlp2*mm";
			}
			elsif($outsheet==2 || $outsheet==3){
				$detector{"dimensions"} = "$out_tile2_pDz*mm $out_tile2_pTheta*deg -$out_tile2_pPhi*deg $out_tile2_pDy1*mm $out_tile2_pDx1*mm $out_tile2_pDx2*mm $out_tile2_pAlp1*deg $out_tile2_pDy2*mm $out_tile2_pDx3*mm $out_tile2_pDx4*mm $out_tile2_pAlp2*mm";
			}
			$detector{"material"} = $out_scint_mat;
			if($outsheet==1 || $outsheet==3){ $detector{"rotation"} = "0*deg 0*deg 180*deg";}
			elsif($outsheet==0 || $outsheet==2){ $detector{"rotation"} = "0*deg 0*deg 0*deg";}
			$detector{"visible"} = "1";
			$detector{"sensitivity"} = "no";
			$detector{"hit_type"}    = "flux";
			$detector{"identifiers"} = "no";
			print_det(\%configuration, \%detector);
		}
	}
}



sub build_proto()
{
	build_box();
	build_hcal_in();
	build_hcal_out();
}
