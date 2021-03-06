$MVC.Native ={};
$MVC.Native.extend = function(class_name, source){
	if(!$MVC[class_name]) $MVC[class_name] = {};
	var dest = $MVC[class_name];
	for (var property in source){
		dest[property] = source[property];
		if(!$MVC._no_conflict){
			window[class_name][property] = source[property];
			if(typeof source[property] == 'function'){
				var names = source[property].toString().match(/^[\s\(]*function[^(]*\((.*?)\)/)[1].split(",");
    			if( names.length == 1 && !names[0]) continue;
				var first_arg = names[0].replace(/^\s+/, '').replace(/\s+$/, '');
				if( first_arg.match(class_name.toLowerCase()) || (first_arg == 'func' && class_name == 'Function' )  ){
					$MVC.Native.set_prototype(class_name, property, source[property]);
				}
			}
		}
	}
};
$MVC.Native.set_prototype = function(class_name, property_name, func){
	window[class_name].prototype[property_name] = function(){
		var args = [this];
		for (var i = 0, length = arguments.length; i < length; i++) args.push(arguments[i]);
		return func.apply(this,args  );
	};
};

//Object helpers
$MVC.Object = {};
$MVC.Object.extend = function(destination, source) {
  for (var property in source)
    destination[property] = source[property];
  return destination;
};
//these are really only for forms
$MVC.Object.to_query_string = function(object,name){
	if(!object) return null;
	if(typeof object == 'string') return object;
	return $MVC.Object.to_query_string.worker(object,name).join('&');
};
$MVC.Object.to_query_string.worker = function(obj,name){
	var parts2 = [];
	for(var thing in obj){
		if(obj.hasOwnProperty(thing)) {
			var value = obj[thing];
			if(typeof value != 'object'){
				var nice_val = encodeURIComponent(value.toString());
				var newer_name = encodeURIComponent(name ? name+'['+thing+']' : thing) ;
				parts2.push( newer_name+'='+nice_val )  ;
			}else{
				if(name){
					parts2 = parts2.concat( $MVC.Object.to_query_string.worker(value, name+'['+thing+']')   );
				}else{
					parts2 = parts2.concat( $MVC.Object.to_query_string.worker(value, thing)  );
				}
			}
		}
	}
	return parts2;
};

//String Helpers
$MVC.Native.extend('String', {
	capitalize : function(string) {
		return string.charAt(0).toUpperCase()+string.substr(1).toLowerCase();
	},
	include : function(string, pattern){
		return string.indexOf(pattern) > -1;
	},
	ends_with : function(string, pattern) {
	    var d = string.length - pattern.length;
	    return d >= 0 && string.lastIndexOf(pattern) === d;
	},
	camelize: function(string){
		var parts = string.split(/_|-/);
		for(var i = 1; i < parts.length; i++)
			parts[i] = $MVC.String.capitalize(parts[i]);
		return parts.join('');
	}
});

//Date Helpers, probably should be moved into its own class

//Array Helpers
$MVC.Native.extend('Array',{ 
	include: function(array, thing){
		for(var i=0; i< array.length; i++){
			if(array[i] == thing) return true;
		}
		return false;
	},
	from: function(iterable){
		 if (!iterable) return [];
		var results = [];
	    for (var i = 0, length = iterable.length; i < length; i++)
	      results.push(iterable[i]);
	    return results;
	}
});
//Function Helpers
$MVC.Native.extend('Function', {
	bind: function(func) {
	  var args = $MVC.Array.from(arguments);
	  args.shift();args.shift();
	  var __method = func, object = arguments[1];
	  return function() {
	    return __method.apply(object, args.concat($MVC.Array.from(arguments) )  );
	  }
	}
});




