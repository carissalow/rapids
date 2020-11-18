window.addEventListener("DOMContentLoaded", function() {
  
    var xhr = new XMLHttpRequest();
    xhr.open("GET", window.location + "../versions.json");
    xhr.onload = function() {
      var versions = JSON.parse(this.responseText);
  
      var realVersion = versions.find(function(i) {
        return i.version === CURRENT_VERSION ||
               i.aliases.includes(CURRENT_VERSION);
      }).version;

      console.log(version)
      console.log(CURRENT_VERSION)
    };
    xhr.send();
  });
  