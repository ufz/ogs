var baseurl = "https://discourse.opengeosys.org";
var announcementsurl = baseurl + "/c/announcements.json";
var usabilityurl = baseurl + "/search.json?q=category%3AUsability";

function newsEntry(element, i, data) {
  var string = "<li>";
  string += "<a href=" + baseurl + "/t/" + element.slug + ">";
  if (i === 0) {
    string += "<i class='fas fa-plus-circle'></i> ";
  }
  else {
    string += "<i class='fal fa-plus-circle'></i> ";
  }
  string += element.title;
  string += "</a>";
  // avatar_template
  // https://discourse.opengeosys.org/user_avatar/discourse.opengeosys.org/bilke/90/405_2.png
  string += element.created_at + " | " + element.tags + " | " + data.users.find(author => author.id === element.posters[0].user_id).name
  string += "</li>";

  return string;
}

if (location.pathname == "/") {
  $.getJSON(announcementsurl, function (data) {
    data.topic_list.topics.slice(0, 3).forEach(function (element, i) {
      $("#news").append(newsEntry(element, i, data));
    });
  });

  $.getJSON(usabilityurl, function (data) {
    data.topics.slice(0, 3).forEach(function (element, i) {
      $("#discussions").append(newsEntry(element, i));
    });
  });
}
