---
layout: post
title: Jekyll+GitHub Pages搭建个人博客
tags: [website, blog]
cover-img: /assets/img/covers/technology01.png
thumbnail-img: /assets/img/mds/octojekyllx.png
music: /assets/audio/friendships.mp3
---

Jekyll是一个静态网站生成器，基于Ruby语言实现。它能够将一定格式的文本文件（如：MarkDown）转换成静态的HTML页面，并提供Liquid模板引擎进行页面渲染；生成的静态网站可以发布到GitHub Pages等网站托管，从而实现自己的项目页面或者个人博客等。此外，Jekyll有丰富的三方开源的个性化、多用途主题。

### Jeykll搭建新站点

以Windows 10为例，搭建Jekyll网站：

1. 安装Jekyll依赖  
下载和安装[Ruby](https://rubyinstaller.org/downloads/)，注意安装时将Ruby添加到系统环境变量。

1. 安装jekyll和bundler gems

```shell
gem install jekyll bundler
```

1. 创建新站点

```shell
jekyll new myblog, 将会在运行目录中生成myblog目录，作为站点根目录  
jekyll new --skip-bundle ., 将以运行目录作为站点根目录
```

新站点目录的内容为：  
**_posts:** 目录，存放写作内容文档，文档格式为日期-名称-后缀，如：2018-08-20-xxxxx.md  
**.gitignore:** 文件，  
**_config.yml:** 文件，站点的配置文件，配置站点的名称、主题风格等  
**404.html:** 文件，404返回页面  
**about.markdown:** 文件，介绍文档  
**Gemfile:**  文件，站点依赖配置文件，配置依赖的库文件等，如使用GitHub Pages托管，需要在此处配置  
**Gemfile.lock:** 文件，  
**index.markdown:** 文件，站点的主页面  
另外，可以可以根据需要添加目录或者页面文档：  
**_layouts:** 目录，存放网页包装模板  
**_includes:** 目录，存放导航栏、页脚等设置文件  
**_data:** 目录，数据文件(支持YAML、JSON、CSV)  
**assets:** 目录，存放静态文件(包括CSS、JS、images等)  

1. 安装站点依赖的库文件，即Gemfile的配置

```shell
bundle add webrick  #Ruby版本3.0.0或以上需要运行此处
bundle install
```

1. 本地站点测试

```shell
bundle exec jekyll serve
```

将会在站点根目录下生成 **_site** 目录，通过http://localhost:4000可访问站点

### 主题应用

Jekyll 有丰富的三方[主题](https://github.com/topics/jekyll-theme)支持，[Beautiful Jekyll](https://github.com/daattali/beautiful-jekyll) 是一个即用型模板，可快速创建漂亮的网站，非常适合个人网站、博客或简单的项目网站，并有详细的使用说明。

beautiful-jekyll写作文档（md
文件）的常用参数有：

|Parameter|Description|
|--|--|
|title|Page or blog post title|
|subtitle|Short description of page or blog post that goes under the title|
|tags|List of tags to categorize the post. Separate the tags with commas and place them inside square brackets. Example: [personal, analysis, finance]|
|cover-img|Include a large full-width image at the top of the page. You can either provide the path to a single image (eg. "/path/to/img") , or a list of images to cycle through (eg. ["/path/img1", "/path/img2"]). If you want to add a caption to an image, then you must use the list notation (use [] even if you have only one image), and each image should be provided as "/path/to/img" : "Caption of image".|
|thumbnail-img|For blog posts, if you want to add a thumbnail that will show up in the feed, use thumbnail-img: /path/to/image. If no thumbnail is provided, then cover-img will be used as the thumbnail. You can use thumbnail-img: "" to disable a thumbnail.|
|music|页面音乐，仅在含cover-image的页面生效，需要提供音乐路径如："/path/to/music"|

### GitHub Pages 展示

GitHub Pages 可以用来展示一些开源项目、主持博客甚或分享个人简历等。完成站点构建后，将站点根目录推送至GitHub仓库，并设置Pages后即能展示。详见[Pages文档](https://docs.github.com/cn/pages/setting-up-a-github-pages-site-with-jekyll/creating-a-github-pages-site-with-jekyll)。要注意的是仓库名需要以.github.io作为后缀。

### 美化

- 添加背景音乐

1、音乐外链，以QQ音乐为例，根据音乐分享的链接，在浏览器中获取歌曲id，如101167817
<iframe frameborder="0" border="0" marginwidth="0" marginheight="0" width=240 height=80
src="//i.y.qq.com/n2/m/outchain/player/index.html?songid=101167817">
</iframe>

代码为：

```html
<iframe frameborder="0" border="0" marginwidth="0" marginheight="0" width=240 height=80
src="//i.y.qq.com/n2/m/outchain/player/index.html?songid=101167817">
</iframe>
```

2、audio元素，结合css控制播放按钮的样式
<audio controls loop>
  <source src="https://zer0o0.github.io/intergap.github.io/assets/audio/silence_voice.mp3" type="audio/mpeg">
  Your browser does not support this audio format.
</audio>

代码为：

```html
<audio controls loop>
  <source src="https://zer0o0.github.io/intergap.github.io/assets/audio/silence_voice.mp3" type="audio/mpeg">
  Your browser does not support this audio format.
</audio>
```

### 参考

[GitHub Pages指南](https://docs.github.com/cn/pages/quickstart)  
[Jekyll官网](https://jekyllrb.com/)  
[beautiful-jekyll theme文档](https://github.com/daattali/beautiful-jekyll)
