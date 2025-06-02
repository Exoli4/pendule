const menus = document.querySelectorAll("div[id^='menu']");
const stateEnum = Object.freeze({
    T2_NORM_SANS: 0,
    T2_EXTREME_SANS: 1,
    T2_2D: 2,
    T2_AVEC: 3
});

let calcState;


window.onload = () =>{
    if(history.state === null)
        showMenu(0);
    else
        showMenu(history.state);
}

window.onpopstate = () =>{
    showMenu(history.back());
}

function showMenu(num){
    menus.forEach(menu =>{
         menu.hidden = true; 
         menu.className = "";
        }
        );

    try{
        menus[num].hidden = false;
        menus[num].className += "flexible";
    }catch (err){
        menus[history.back()].hidden = false;
        menus[history.back()].classList += "flexible";
    }

    if(num === 3){
            let titre = document.querySelector("#title");
            switch(calcState){
                case 0: titre.textContent = "Comparaison T2 (pas de force externe)";
                    break;
                case 1: titre.textContent = "Comparaison T2 sans extrémités (pas de force externe)";
                    break;
                case 2: titre.textContent = "Comparaison T2 2D (pas de force externe)";
                    break;
                case 3: titre.textContent = "Comparaison T2 (pas de force externe)";
                        const form = document.querySelector("form");
                        const div = document.createElement("div");

                        div.innerHTML +=("<label for=\"force\">Force (N): </label>");
                        div.innerHTML +=("<input type=\"number\" placeholder=\"1.0\" step=\"0.1\" id=\"force\">");

                        form.insertBefore(div, form.lastElementChild);
                    break;
                default: titre.textContent = "null";
        }
    }

    history.pushState(num, "", "");
}

function setCalcState(num){
    switch(num){
        case 0: calcState = stateEnum.T2_NORM_SANS;
            break;
        case 1: calcState = stateEnum.T2_EXTREME_SANS;
            break;
        case 2: calcState = stateEnum.T2_2D;
            break;
        case 3: calcState = stateEnum.T2_AVEC;
            break;
        default: console.log("Vaue not valid");
        break;
    }

    showMenu(3);
}