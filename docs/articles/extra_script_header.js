if (localStorage.getItem('passwordProvided') === null){
    function run(){
        var password = prompt("Enter Password");
        while (password != 'metapsyTools123!'){
            alert('Wrong Password');
            var password = prompt("Enter Password");
        }
        localStorage.setItem('passwordProvided', 'true');
    }
    run();
}
