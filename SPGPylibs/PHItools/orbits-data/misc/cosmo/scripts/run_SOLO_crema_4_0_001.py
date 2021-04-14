 #
# This example script makes the program take users for a 2 minute trip
# to the Earth and Moon.
#
# Import Cosmographia scripting module.
#
import cosmoscripting
cosmo = cosmoscripting.Cosmo()

waittime    = 1
msgwaittime = 1
msgshowtime = 2

testjson = "/Users/mcosta/SOLAR-ORBITER/misc/cosmo/scenarios/load_SOLO_default_001.json"

# cosmo.displayNote('cosmo.loadCatalogFile(testjson)', msgshowtime).wait(0)
cosmo.loadCatalogFile(testjson).wait( waittime )

#
# Set time using setTime to 2015-10-31 23:50:00 UTC. The total time to
# run the commands below will be 3 seconds -- 2 in wait + another 1 in
# wait. The setTime command is instantaneous.
#
# cosmo.displayNote( "Set time to 2025-10-27 23:50 UTC", 2 ).wait( 0 )
cosmo.setTime( "2020-04-01 04:30:00.000 UTC" )

cosmo.gotoObject( "SOLO", 2 )

#cosmo.rollLeft(90.0, 3.0).wait( waittime )

cosmo.faster10x()
cosmo.faster10x()

#
# Reset all visual attributes and settings that we changed.
# All commands are instanteneous and will occur during the 3
# second wait set by wait.
#
cosmo.showToolBar()
cosmo.showStatusMessages()
cosmo.showInfoText()
cosmo.showNormalWindow()


#
# End of the example script.
#
