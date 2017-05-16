#import subprocess # just to call an arbitrary command e.g. 'ls'
#
## enter the directory like this:
#with cd("C:\Users\Zahra\Documents\GitHub\SurfaceRendering\mitsuba\dist"):
#   # we are in ~/Library
#   subprocess.call("ls")



# Import system modules
#import os
#os.chdir("C:\Users\Zahra\Documents\GitHub\SurfaceRendering\mitsuba\dist");
#print(os.getcwd())
#os.system('mitsuba C:\Users\Zahra\Documents\GitHub\SurfaceRendering\SurfaceRendering\Output\diff_angle_height_flyaways\gabardine_diff_config_height_flyaways.xml')
##os.system("mitsuba")

import subprocess
import os
#os.chdir("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/mitsuba/dist");
#subprocess.call('mitsuba -D x="0,7,9" C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/diff_angle_height_flyaways/gabardine_diff_config_height_flyaways.xml');
#os.chdir("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/diff_angle_height_flyaways");
#os.rename("gabardine_diff_config_height_flyaways.exr", "origin_079.exr")              


os.chdir("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/mitsuba/dist");
subprocess.call('mitsuba -D x="0,10,7" C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/diff_angle_height_flyaways/gabardine_diff_config_height_flyaways.xml');
os.chdir("C:/Users/Zahra/Documents/GitHub/SurfaceRendering/SurfaceRendering/Output/diff_angle_height_flyaways");
os.rename("gabardine_diff_config_height_flyaways.exr", "origin_0107.exr")  