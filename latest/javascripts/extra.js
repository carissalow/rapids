window.addEventListener("DOMContentLoaded", function() {
    if(!window.location.pathname.endsWith("/latest/")){
        document.querySelector("div[data-md-component=announce]").innerHTML = "<div id='announce-msg'>You are seeing the docs for a previous version of RAPIDS, <a href='latest/'>click here to go to the latest</a></div>"
    }
  });
