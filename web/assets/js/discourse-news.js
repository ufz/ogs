var baseurl = "https://discourse.opengeosys.org";
var announcementsurl = baseurl + "/search.json?q=category%3AAnnouncements";
var usabilityurl = baseurl + "/search.json?q=category%3AUsability";

function newsEntry(element, i) {
  var string = "<li>";
  string += "<a href=" + baseurl + "/t/" + element.slug + ">";
  if(i === 0) {
    string += "<i class='fas fa-plus-circle'></i> ";
  }
  else {
    string += "<i class='fal fa-plus-circle'></i> ";
  }
  string += element.title;
  string += "</a>";
  string += "</li>";

  return string;
}

$.getJSON(announcementsurl, function(data) {
  data.topics.slice(0,3).forEach(function (element, i) {
    $("#news").append(newsEntry(element, i));
  });
});

$.getJSON(usabilityurl, function(data) {
  data.topics.slice(0,3).forEach(function (element, i) {
    $("#discussions").append(newsEntry(element, i));
  });
});
