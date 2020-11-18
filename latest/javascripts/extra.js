xhr.onload = function() { 
    var versions = JSON.parse(this.responseText); 
 
    var realVersion = versions.find(function(i) { 
      return i.version === CURRENT_VERSION || 
             i.aliases.includes(CURRENT_VERSION); 
    }).version; 
    console.log(realVersion)
    console.log(CURRENT_VERSION)
}

console.log("hi")