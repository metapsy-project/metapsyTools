const footerRight = document.getElementsByClassName('pkgdown-footer-right')[0];

var imgTUM = document.createElement('img');
imgTUM.src = 'assets/figures/tum.png';
imgTUM.width = "100";
imgTUM.style = "Padding: 20px 20px 20px 20px;";

var imgVU = document.createElement('img');
imgVU.src = 'assets/figures/vu.png';
imgVU.width = "125";
imgVU.style = "Padding: 20px 20px 20px 20px;";

footerRight.appendChild(imgTUM);
footerRight.appendChild(imgVU);

const navbarRight = document.getElementsByClassName('navbar-nav')[1];

var imgTumWhiteContainer = document.createElement('li');
var imgTumWhite = document.createElement('img');
imgTumWhite.src = 'figures/white_tum.png';
imgTumWhite.height = "30";
imgTumWhite.style = "Padding: 10px 0px 0px 10px;";
imgTumWhiteContainer.appendChild(imgTumWhite);
navbarRight.appendChild(imgTumWhiteContainer);

var imgVuWhiteContainer = document.createElement('li');
var imgVuWhite = document.createElement('img');
imgVuWhite.src = 'figures/white_vu.png';
imgVuWhite.height = "30";
imgVuWhite.style = "Padding: 10px 0px 0px 15px;";
imgVuWhiteContainer.appendChild(imgVuWhite);
navbarRight.appendChild(imgVuWhiteContainer);
