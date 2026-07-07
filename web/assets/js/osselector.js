// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause
jQuery.fn.visible = function () {
  return this.css('visibility', 'visible');
};

jQuery.fn.invisible = function () {
  return this.css('display', 'none');
};

$("#btn-win").click(function () { selectOS("win"); });
$("#btn-linux").click(function () { selectOS("linux"); });
$("#btn-mac").click(function () { selectOS("mac"); });

function selectOS(os) {
  var os_list = ['win', 'linux', 'mac'];
  for (var i = 0; i < os_list.length; i++) {
    var current_os = os_list[i];
    if (current_os == os) {
      $("." + current_os).show();
    }
    else {
      $("." + current_os).hide();
    }
    $("#btn-" + current_os).removeClass("active");
  }
  $("[data-os]").hide();
  $("[data-os~='all'], [data-os~='" + os + "']").show();
  $("#btn-" + os).addClass("active");
  window.localStorage.setItem("selectedOS", os);
}

function hasOS(os) {
  return $("." + os).length > 0 || $("[data-os~='all']").length > 0 || $("[data-os~='" + os + "']").length > 0;
}

var available_os = ['win', 'linux', 'mac'].filter(function (os) { return hasOS(os); });

for (var i = 0; i < ['win', 'linux', 'mac'].length; i++) {
  var os = ['win', 'linux', 'mac'][i];
  if (hasOS(os)) {
    $("#btn-" + os).show();
  }
  else {
    $("#btn-" + os).hide();
  }
}

if (available_os.length > 0) {
  $("#os-selector").visible();
  var os = window.localStorage.getItem("selectedOS");
  if (os && available_os.indexOf(os) !== -1) {
    $("#btn-" + os).click();
  }
  else {
    $("#btn-" + available_os[0]).click();
  }
}
else {
  $("#os-selector").invisible();
}
