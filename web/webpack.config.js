var webpack = require('webpack');
const ExtractTextPlugin = require('extract-text-webpack-plugin')
const path = require('path');

const sharp = require('sharp');

module.exports = {
  mode: 'development',
  entry: './src/index.js',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'bundle.js',
  },
  module: {
    rules: [
      {
        test: /\.css$/,
        use: ExtractTextPlugin.extract({
          fallback: 'style-loader',
          use: [
            { loader: 'css-loader', options: { importLoaders: 1 } },
            'postcss-loader'
          ]
        })
      },
      {
        test: /\.(gif|jpe?g|png|tiff)(\?.*)?$/,
        use: [
          {
            loader: 'sharp-loader',
            query: {
              context: path.join(__dirname, 'content'),
              name: '[path][name].[ext].[preset]',
              cacheDirectory: true,
              presets: {
                thumbnail: {
                  width: 320,
                  height: 200,
                  quality: 70,
                }
              }
            },
          }
        ]
      },
      {
        test: /\.(eot|svg|ttf|woff|woff2)$/,
        use: [
          {
            loader: 'file-loader',
            options: {}
          }
        ]
      }
    ]
  },
  plugins: [
    new ExtractTextPlugin('styles.css'),
    new webpack.ProvidePlugin({
      $: "jquery",
      jQuery: "jquery",
      'window.jQuery': 'jquery'
    })
  ],
  devtool: 'source-map'
}
