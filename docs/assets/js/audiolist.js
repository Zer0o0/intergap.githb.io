window.onload = function(){
    var musicList = ["/assets/audio/silence_voice.mp3",]  //音乐列表
    playMusic(musicList);
  }

  function playMusic(musicList){ 
    var player = new Audio(); 
    player.preload = true;
    player.controls = true; 
    var src = musicList.pop();
    player.src = src;
    musicList.unshift(src);
    player.addEventListener("ended",playEndedHandler,false);
    player.play();
    document.getElementById("audioBox").appendChild(player);
    player.loop = false

    function playEndedHandler(){ 
      src = musicList.pop(); 
      player.src = src; 
      musicList.unshift(src); 
      player.play();
    } 
  } 

//HTML 文档中使用js
//<div id="audioBox"> 
//  <script type="text/javascript" src="/assets/js/audiolist.js"> </script> 
//</div>
