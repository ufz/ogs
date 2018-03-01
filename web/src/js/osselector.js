jQuery.fn.visible = function() {
    return this.css('visibility', 'visible');
};

jQuery.fn.invisible = function() {
    return this.css('visibility', 'hidden');
};

$("#btn-win").click(function(){ selectOS("win"); });
$("#btn-linux").click(function(){ selectOS("linux"); });
$("#btn-mac").click(function(){ selectOS("mac"); });

function selectOS(os) {
  var os_list = ['win', 'linux', 'mac'];
  for (var i = 0; i < os_list.length; i++) {
    var current_os = os_list[i];
    if (current_os == os) {
      $("." + current_os).show();
      $("#btn-" + current_os).addClass("active");
    }
    else {
      $("." + current_os).hide();
      $("#btn-" + current_os).removeClass("active");
    }
  }
  window.localStorage.setItem("selectedOS", os);
}

if ($(".win").length > 0) {
  $("#os-selector").visible();
  var os = window.localStorage.getItem("selectedOS");
  if (os) {
    $("#btn-" + os).click();
  }
  else {
    $("#btn-win").click();
  }
}
