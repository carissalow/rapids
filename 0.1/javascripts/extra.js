window.addEventListener("DOMContentLoaded", function() {
    var xhr = new XMLHttpRequest();
    xhr.open("GET", window.location.origin + "/versions.json");
    xhr.onload = function() {
      var versions = JSON.parse(this.responseText);
      latest_version = ""
      for(id in versions)
          if(versions[id]["aliases"].length > 0 && versions[id]["aliases"].includes("latest")) 
              latest_version = "/" + versions[id].version + "/"
      if(!window.location.pathname.includes("/latest/") && (latest_version.length > 0 && !window.location.pathname.includes(latest_version)))
          document.querySelector("div[data-md-component=announce]").innerHTML = "<div id='announce-msg'>You are seeing the docs for a previous version of RAPIDS, <a href='http://www.rapids.science/'>click here to go to the latest</a></div>"
    };
    xhr.send();
  });
