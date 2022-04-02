const footerRight = document.getElementsByClassName('pkgdown-footer-right')[0];

var imgTUM = document.createElement('img');
imgTUM.src = '../assets/figures/tum.png';
imgTUM.width = "100";
imgTUM.style = "Padding: 20px 20px 20px 20px;";

var imgVU = document.createElement('img');
imgVU.src = '../assets/figures/vu.png';
imgVU.width = "125";
imgVU.style = "Padding: 20px 20px 20px 20px;";

footerRight.appendChild(imgTUM);
footerRight.appendChild(imgVU);
