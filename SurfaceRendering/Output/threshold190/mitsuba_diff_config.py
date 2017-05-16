import subprocess
import os

#os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
#subprocess.call('mitsuba -D x="0,7,7" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
#os.rename("untiled_hybrid.exr", "origin_0_7_7.exr")  
#
#os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
#subprocess.call('mitsuba -D x="0,9,7" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
#os.rename("untiled_hybrid.exr", "origin_0_9_7.exr")  
#
#os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
#subprocess.call('mitsuba -D x="0,10,7" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
#os.rename("untiled_hybrid.exr", "origin_0_10_7.exr")  
#
#os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
#subprocess.call('mitsuba -D x="0,7,9" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
#os.rename("untiled_hybrid.exr", "origin_0_7_9.exr") 

## now comment the shape for the volume and run again
#os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
#subprocess.call('mitsuba -D x="0,7,7" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
#os.rename("untiled_hybrid.exr", "height_origin_0_7_7.exr")  
#
#os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
#subprocess.call('mitsuba -D x="0,9,7" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
#os.rename("untiled_hybrid.exr", "height_origin_0_9_7.exr")  
#
#os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
#subprocess.call('mitsuba -D x="0,10,7" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
#os.rename("untiled_hybrid.exr", "height_origin_0_10_7.exr")  
#
#os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
#subprocess.call('mitsuba -D x="0,7,9" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
#os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
#os.rename("untiled_hybrid.exr", "height_origin_0_7_9.exr") 

# now rotate around x axis
os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
subprocess.call('mitsuba -D x="0" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
os.rename("untiled_hybrid.exr", "height_angle_0.exr")  

os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
subprocess.call('mitsuba -D x="10" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
os.rename("untiled_hybrid.exr", "height_angle_10.exr")   

os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
subprocess.call('mitsuba -D x="20" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
os.rename("untiled_hybrid.exr", "height_angle_20.exr")  

os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
subprocess.call('mitsuba -D x="30" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
os.rename("untiled_hybrid.exr", "height_angle_30.exr") 

os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
subprocess.call('mitsuba -D x="40" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
os.rename("untiled_hybrid.exr", "height_angle_40.exr") 

os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
subprocess.call('mitsuba -D x="45" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
os.rename("untiled_hybrid.exr", "height_angle_45.exr") 

os.chdir("E:/zahra/SurfaceRendering-master/mitsuba/dist");
subprocess.call('mitsuba -D x="60" E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190/untiled_hybrid.xml');
os.chdir("E:/zahra/SurfaceRendering-master/SurfaceRendering/Output/threshold190");
os.rename("untiled_hybrid.exr", "height_angle_60.exr") 