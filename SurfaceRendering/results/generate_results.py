from shutil import copyfile
import subprocess
import os

os.chdir("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/mitsuba/dist");

#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/hybrid_vs_ref/170_hybrid_angle_45.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold175/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold175/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/hybrid_vs_ref/175_hybrid_angle_45.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold180/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold180/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/hybrid_vs_ref/180_hybrid_angle_45.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold185/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold185/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/hybrid_vs_ref/185_hybrid_angle_45.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold190/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/hybrid_vs_ref/190_hybrid_angle_45.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold232/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold232/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/hybrid_vs_ref/232_hybrid_angle_45.exr")



##### reference volume 
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/hybrid_vs_ref/reference_volume.exr")



#### cropped volume 
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/170_cropped_vol.exr")

#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold175/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold175/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/175_cropped_vol.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold180/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold180/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/180_cropped_vol.exr")

#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold185/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold185/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/185_cropped_vol.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold190/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/190_cropped_vol.exr")



### heightmap only
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/170_height.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold175/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold175/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/175_height.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold180/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold180/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/180_height.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold185/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold185/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/185_height.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold190/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/190_height.exr")
#
#subprocess.call('mitsuba C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold232/untiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold232/untiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/232_height.exr")



##### zoomed-in results from tiled hybrid
#subprocess.call('mitsuba -D tileX="6" -D tileY="6" -D scaleX="1" -D scaleY="1" -D fov="5" C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/tiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/tiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/zoomed_tiled/6x6tiled_hybrid_170.exr");
#        
subprocess.call('mitsuba -D tileX="6" -D tileY="6" -D scaleX="1" -D scaleY="1" -D fov="5" C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold175/tiled_hybrid.xml');
copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold175/tiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/zoomed_tiled/6x6tiled_hybrid_175.exr");

subprocess.call('mitsuba -D tileX="6" -D tileY="6" -D scaleX="1" -D scaleY="1" -D fov="5" C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold180/tiled_hybrid.xml');
copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold180/tiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/zoomed_tiled/6x6tiled_hybrid_180.exr");

subprocess.call('mitsuba -D tileX="6" -D tileY="6" -D scaleX="1" -D scaleY="1" -D fov="5" C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold185/tiled_hybrid.xml');
copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold185/tiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/zoomed_tiled/6x6tiled_hybrid_185.exr");

subprocess.call('mitsuba -D tileX="6" -D tileY="6" -D scaleX="1" -D scaleY="1" -D fov="5" C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold190/tiled_hybrid.xml');
copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold190/tiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/zoomed_tiled/6x6tiled_hybrid_190.exr");

subprocess.call('mitsuba -D tileX="6" -D tileY="6" -D scaleX="1" -D scaleY="1" -D fov="5" C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold232/tiled_hybrid.xml');
copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold232/tiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/zoomed_tiled/6x6tiled_hybrid_232.exr");        
#
#subprocess.call('mitsuba -D tileX="10" -D tileY="15" -D scaleX="1.67" -D scaleY="2.5" -D fov="2" C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/tiled_hybrid.xml');
#copyfile("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/threshold170/tiled_hybrid.exr", "C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/results/height_vs_cropVol/test.exr")
