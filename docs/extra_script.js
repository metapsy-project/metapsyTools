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

var nav = document.getElementsByTagName('nav');
var strip = document.createElement('div');
strip.id = "strip";
strip.style.backgroundColor = "red";
strip.style.position = "absolute";
strip.style.width = "100%";
strip.style.marginTop = "-50px";
strip.style.fontSize = "13px";
strip.style.paddingBottom = "10px";
strip.style.paddingLeft = "28px";
strip.style.backgroundImage = "linear-gradient(90deg, rgba(94, 186, 189, 1), rgba(3, 114, 183, 1))";
strip.style.marginTop = "-78px";
var container = document.getElementsByClassName('container');
container[0].style.paddingTop = "13px";
nav[0].insertBefore(strip, nav[0].firstChild);

