/**
 * Event describes 2 functions
 * 	Event.observe
 * 	Event.stopObserving
 */

if(typeof Prototype == 'undefined') 
	include("standard");
else
	$MVC.Event = Event;
//jQuery's wont work for controllers because it doesn't allow capture