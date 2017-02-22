// load all plugins in 'devDependencies' into the variable $
var gulp = require('gulp');

const $ = require('gulp-load-plugins')({
        rename: {
            'gulp-add-src': 'add_src',
            'merge-stream': 'merge_stream',
            'gulp-clean-css': 'clean_css',
            'gulp-uglifycss': 'uglifycss'
        },
        pattern: ['*'],
        scope: ['devDependencies']
    });

// package vars
const pkg = require('./package.json');

// scss - build the scss to the build folder, including the required paths, and writing out a sourcemap
gulp.task('scss', () => {
    $.fancyLog("-> Compiling scss: " + pkg.paths.dist.css + pkg.vars.scssName);
    sassStream = gulp.src(pkg.paths.src.scss + pkg.vars.scssName)
        // .pipe($.plumber({ errorHandler: onError }))
        // .pipe($.sourcemaps.init())
        .pipe($.sass({
                includePaths: [
                  pkg.paths.scss,
                  "./node_modules"
                ]
            })
            .on('error', $.sass.logError))
        .pipe($.cached('sass_compile'))
        .pipe($.autoprefixer())
        // .pipe($.sourcemaps.write('./'))
        .pipe($.size({ showFiles: true }));

    cssStream = gulp.src(pkg.globs.distCss);

    return $.merge_stream(sassStream, cssStream)
        .pipe($.concat(pkg.vars.cssName))
        .pipe($.clean_css({compatibility: 'ie8', level: 1}))
        .pipe($.size({title: 'Cleaned', showFiles: true}))
        .pipe(gulp.dest(pkg.paths.dist.css));
});


gulp.task('watch', function() {
  gulp.watch(pkg.paths.src.scss + pkg.vars.scssName, ['scss']);
  gulp.watch('./package.json', ['scss']);
});

gulp.task('clean', function() {
    return $.del([
        'static',
        'public',
        'content/internal/news.md',
        'data/news.json'
    ]);
});

gulp.task('clean-all', function() {
    return $.del([
        'static',
        'public',
        'node_modules',
        'content/internal/news.md',
        'data/news.json'
    ]);
});

gulp.task('build', ['scss'])

gulp.task('default', ['scss', 'watch'])
